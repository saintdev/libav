/*
 * AAC coefficients encoder
 * Copyright (C) 2008-2009 Konstantin Shishkov
 *
 * This file is part of Libav.
 *
 * Libav is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * Libav is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Libav; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

/**
 * @file
 * AAC coefficients encoder
 */

/***********************************
 *              TODOs:
 * speedup quantizer selection
 * add sane pulse detection
 ***********************************/

#include "libavutil/libm.h" // brought forward to work around cygwin header breakage

#include <float.h>
#include "libavutil/mathematics.h"
#include "avcodec.h"
#include "put_bits.h"
#include "aac.h"
#include "aacenc.h"
#include "aactab.h"

/** bits needed to code codebook run value for long windows */
static const uint8_t run_value_bits_long[64] = {
     5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,
     5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5, 10,
    10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
    10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 15
};

/** bits needed to code codebook run value for short windows */
static const uint8_t run_value_bits_short[16] = {
    3, 3, 3, 3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 6, 6, 9
};

static const uint8_t *run_value_bits[2] = {
    run_value_bits_long, run_value_bits_short
};


/**
 * Quantize one coefficient.
 * @return absolute value of the quantized coefficient
 * @see 3GPP TS26.403 5.6.2 "Scalefactor determination"
 */
static av_always_inline int quant(float coef, const float Q)
{
    float a = coef * Q;
    return sqrtf(a * sqrtf(a)) + 0.4054;
}

static void quantize_bands(int *out, const float *in, const float *scaled,
                           int size, float Q34)
{
    int i;
    double qc;
    for (i = 0; i < size; i++) {
        qc = scaled[i] * Q34;
        out[i] = (int)FFMIN(qc + 0.4054, 8191.0);
        if (in[i] < 0.0f)
            out[i] = -out[i];
    }
}

static void invquant_bands(float *out, const int *in, int size, float IQ)
{
    int i;

    for (i = 0; i < size; i++) {
        out[i] = (float)((double)in[i] * cbrt(in[i]) * (double)IQ);
        if (in[i] < 0)
            out[i] = -out[i];
    }
}

void ff_aac_quantize_band(int *out, const float *in, const float *scaled, int sf_idx, int size)
{
    const float Q34 = ff_aac_pow34sf_tab[POW_SF2_ZERO - sf_idx + SCALE_ONE_POS - SCALE_DIV_512];

    quantize_bands(out, in, scaled, size, Q34);
}

static void abs_pow34_v(float *out, const float *in, const int size)
{
#ifndef USE_REALLY_FULL_SEARCH
    int i;
    for (i = 0; i < size; i++) {
        float a = fabsf(in[i]);
        out[i] = sqrtf(a * sqrtf(a));
    }
#endif /* USE_REALLY_FULL_SEARCH */
}

void ff_abs_pow34_v(float *out, const float *in)
{
    abs_pow34_v(out, in, 1024);
}

static const uint8_t aac_cb_range [12] = {0, 3, 3, 3, 3, 9, 9, 8, 8, 13, 13, 17};
static const uint8_t aac_cb_maxval[12] = {0, 1, 1, 2, 2, 4, 4, 7, 7, 12, 12, 16};

#define SFB_FUNC_TEMPLATE(type, cb, BT_UNSIGNED, BT_PAIR, BT_ESC) \
static void sfb_ ##type ##_hcb_ ##cb(AACEncContext *s, const int *quants, int size) \
{ \
    const int dim = BT_PAIR ? 2 : 4; \
    const int off = BT_UNSIGNED ? 0 : aac_cb_maxval[cb]; \
\
    if (!cb) \
        return; \
\
    for (int i = 0; i < size; i += dim) { \
        int curidx = 0; \
\
        for (int j = 0; j < dim; j++) { \
            int t; \
            if (BT_UNSIGNED) \
                t = FFMIN(abs(quants[i + j]), aac_cb_maxval[cb]); \
            else \
                t = av_clip(quants[i + j], -aac_cb_maxval[cb], aac_cb_maxval[cb]); \
            curidx *= aac_cb_range[cb]; \
            curidx += t + off; \
        } \
        put_bits_template(s, ff_aac_spectral_bits[cb - 1][curidx], ff_aac_spectral_codes[cb - 1][curidx]); \
        if (BT_UNSIGNED) \
            for (int j = 0; j < dim; j++) \
                if (ff_aac_codebook_vectors[cb - 1][curidx * dim + j] != 0.0f) \
                    put_bits_template(s, 1, quants[i + j] < 0); \
        if (BT_ESC) { \
            for (int j = 0; j < 2; j++) { \
                if (ff_aac_codebook_vectors[cb - 1][curidx * 2 + j] == 64.0f) { \
                    int coef = abs(quants[i + j]); \
                    int len = av_log2(coef); \
\
                    put_bits_template(s, len - 4 + 1, (1 << (len - 4 + 1)) - 2); \
                    put_bits_template(s, len, coef & ((1 << len) - 1)); \
                } \
            } \
        } \
    } \
}

#define SFB_ENCODE_FUNC(cb, BT_UNSIGNED, BT_PAIR, BT_ESC) \
    SFB_FUNC_TEMPLATE(encode, cb, BT_UNSIGNED, BT_PAIR, BT_ESC)

#define SFB_SIZE_FUNC(cb, BT_UNSIGNED, BT_PAIR, BT_ESC) \
    SFB_FUNC_TEMPLATE(size, cb, BT_UNSIGNED, BT_PAIR, BT_ESC)

#define put_bits_template(s, count, value) put_bits(&s->pb, count, value)
SFB_ENCODE_FUNC( 0, 0, 0, 0 )
SFB_ENCODE_FUNC( 1, 0, 0, 0 )
SFB_ENCODE_FUNC( 2, 0, 0, 0 )
SFB_ENCODE_FUNC( 3, 1, 0, 0 )
SFB_ENCODE_FUNC( 4, 1, 0, 0 )
SFB_ENCODE_FUNC( 5, 0, 1, 0 )
SFB_ENCODE_FUNC( 6, 0, 1, 0 )
SFB_ENCODE_FUNC( 7, 1, 1, 0 )
SFB_ENCODE_FUNC( 8, 1, 1, 0 )
SFB_ENCODE_FUNC( 9, 1, 1, 0 )
SFB_ENCODE_FUNC(10, 1, 1, 0 )
SFB_ENCODE_FUNC(11, 1, 1, 1 )
#undef put_bits_template

#define put_bits_template(s, count, value) do { s->bits_encoded += count; } while (0)
SFB_SIZE_FUNC( 0, 0, 0, 0 )
SFB_SIZE_FUNC( 1, 0, 0, 0 )
SFB_SIZE_FUNC( 2, 0, 0, 0 )
SFB_SIZE_FUNC( 3, 1, 0, 0 )
SFB_SIZE_FUNC( 4, 1, 0, 0 )
SFB_SIZE_FUNC( 5, 0, 1, 0 )
SFB_SIZE_FUNC( 6, 0, 1, 0 )
SFB_SIZE_FUNC( 7, 1, 1, 0 )
SFB_SIZE_FUNC( 8, 1, 1, 0 )
SFB_SIZE_FUNC( 9, 1, 1, 0 )
SFB_SIZE_FUNC(10, 1, 1, 0 )
SFB_SIZE_FUNC(11, 1, 1, 1 )
#undef put_bits_template

static void (*sfb_encode_hcb_arr[])(AACEncContext *s, const int *quants, int size) = {
    sfb_encode_hcb_0, sfb_encode_hcb_1, sfb_encode_hcb_2,  sfb_encode_hcb_3,
    sfb_encode_hcb_4, sfb_encode_hcb_5, sfb_encode_hcb_6,  sfb_encode_hcb_7,
    sfb_encode_hcb_8, sfb_encode_hcb_9, sfb_encode_hcb_10, sfb_encode_hcb_11
};

static void (*sfb_size_hcb_arr[])(AACEncContext *s, const int *quants, int size) = {
    sfb_size_hcb_0, sfb_size_hcb_1, sfb_size_hcb_2,  sfb_size_hcb_3,
    sfb_size_hcb_4, sfb_size_hcb_5, sfb_size_hcb_6,  sfb_size_hcb_7,
    sfb_size_hcb_8, sfb_size_hcb_9, sfb_size_hcb_10, sfb_size_hcb_11
};

void ff_aac_sfb_encode_hcb(AACEncContext *s, const int *quants, int size, int cb)
{
    sfb_encode_hcb_arr[cb](s, quants, size);
}

static int sfb_size_hcb(AACEncContext *s, const int *quants, int size, int cb)
{
    s->bits_encoded = 0;
    sfb_size_hcb_arr[cb](s, quants, size);

    return s->bits_encoded;
}

static float band_ssd(const float *coeffs, const float *iquants, int size)
{
    int i;
    float diff, ssd = 0.0f;

    for (i = 0; i < size; i++) {
        diff = coeffs[i] - iquants[i];
        ssd += diff * diff;
    }

    return ssd;
}

/**
 * Calculate rate distortion cost for quantizing with given codebook
 *
 * @return quantization distortion
 */
static float sfb_rd_cost(struct AACEncContext *s, const float *in, const float *scaled,
                          int size, int scale_idx, int cb, const float lambda)
{
    const float Q34 = ff_aac_pow34sf_tab[POW_SF2_ZERO - scale_idx + SCALE_ONE_POS - SCALE_DIV_512];
    const float IQ  = ff_aac_pow2sf_tab [POW_SF2_ZERO + scale_idx - SCALE_ONE_POS + SCALE_DIV_512];
    float cost = 0;
    float ssd;

    if (!cb) {
        int i;

        for (i = 0; i < size; i++)
            cost += in[i]*in[i];

        return cost * lambda;
    }
    quantize_bands(s->qcoefs, in, scaled, size, Q34);
    invquant_bands(s->iquants, s->qcoefs, size, IQ);
    ssd = band_ssd(in, s->iquants, size);

    cost = ssd * lambda + sfb_size_hcb(s, s->qcoefs, size, cb);

    return cost;
}

static float find_max_val(int group_len, int swb_size, const float *scaled) {
    float maxval = 0.0f;
    int w2, i;
    for (w2 = 0; w2 < group_len; w2++) {
        for (i = 0; i < swb_size; i++) {
            maxval = FFMAX(maxval, scaled[w2*128+i]);
        }
    }
    return maxval;
}

static int find_min_book(float maxval, int sf) {
    float Q = ff_aac_pow2sf_tab[POW_SF2_ZERO - sf + SCALE_ONE_POS - SCALE_DIV_512];
    float Q34 = sqrtf(Q * sqrtf(Q));
    int qmaxval, cb;
    qmaxval = maxval * Q34 + 0.4054f;
    if      (qmaxval ==  0) cb = 0;
    else if (qmaxval ==  1) cb = 1;
    else if (qmaxval ==  2) cb = 3;
    else if (qmaxval <=  4) cb = 5;
    else if (qmaxval <=  7) cb = 7;
    else if (qmaxval <= 12) cb = 9;
    else                    cb = 11;
    return cb;
}

/**
 * structure used in optimal codebook search
 */
typedef struct BandCodingPath {
    int prev_idx; ///< pointer to the previous path point
    float cost;   ///< path cost
    int run;
} BandCodingPath;

/**
 * Encode band info for single window group bands.
 */
static void encode_window_bands_info(AACEncContext *s, SingleChannelElement *sce,
                                     int win, int group_len, const float lambda)
{
    BandCodingPath path[120][12];
    int w, swb, cb, start, size;
    int i, j;
    const int max_sfb  = sce->ics.max_sfb;
    const int run_bits = sce->ics.num_windows == 1 ? 5 : 3;
    const int run_esc  = (1 << run_bits) - 1;
    int idx, ppos, count;
    int stackrun[120], stackcb[120], stack_len;
    float next_minrd = INFINITY;
    int next_mincb = 0;

    start = win*128;
    for (cb = 0; cb < 12; cb++) {
        path[0][cb].cost     = 0.0f;
        path[0][cb].prev_idx = -1;
        path[0][cb].run      = 0;
    }
    for (swb = 0; swb < max_sfb; swb++) {
        size = sce->ics.swb_sizes[swb];
        if (sce->zeroes[win*16 + swb]) {
            for (cb = 0; cb < 12; cb++) {
                path[swb+1][cb].prev_idx = cb;
                path[swb+1][cb].cost     = path[swb][cb].cost;
                path[swb+1][cb].run      = path[swb][cb].run + 1;
            }
        } else {
            float minrd = next_minrd;
            int mincb = next_mincb;
            next_minrd = INFINITY;
            next_mincb = 0;
            for (cb = 0; cb < 12; cb++) {
                float cost_stay_here, cost_get_here;
                float rd = 0.0f;
                for (w = 0; w < group_len; w++) {
                    FFPsyBand *band = &s->psy.ch[s->cur_channel].psy_bands[(win+w)*16+swb];
                    rd += sfb_rd_cost(s, sce->coeffs + start + w*128,
                                       sce->scoeffs + start + w*128, size,
                                       sce->sf_idx[(win+w)*16+swb], cb,
                                       lambda / band->threshold);
                }
                cost_stay_here = path[swb][cb].cost + rd;
                cost_get_here  = minrd              + rd + run_bits + 4;
                if (   run_value_bits[sce->ics.num_windows == 8][path[swb][cb].run]
                    != run_value_bits[sce->ics.num_windows == 8][path[swb][cb].run+1])
                    cost_stay_here += run_bits;
                if (cost_get_here < cost_stay_here) {
                    path[swb+1][cb].prev_idx = mincb;
                    path[swb+1][cb].cost     = cost_get_here;
                    path[swb+1][cb].run      = 1;
                } else {
                    path[swb+1][cb].prev_idx = cb;
                    path[swb+1][cb].cost     = cost_stay_here;
                    path[swb+1][cb].run      = path[swb][cb].run + 1;
                }
                if (path[swb+1][cb].cost < next_minrd) {
                    next_minrd = path[swb+1][cb].cost;
                    next_mincb = cb;
                }
            }
        }
        start += sce->ics.swb_sizes[swb];
    }

    //convert resulting path from backward-linked list
    stack_len = 0;
    idx       = 0;
    for (cb = 1; cb < 12; cb++)
        if (path[max_sfb][cb].cost < path[max_sfb][idx].cost)
            idx = cb;
    ppos = max_sfb;
    while (ppos > 0) {
        cb = idx;
        stackrun[stack_len] = path[ppos][cb].run;
        stackcb [stack_len] = cb;
        idx = path[ppos-path[ppos][cb].run+1][cb].prev_idx;
        ppos -= path[ppos][cb].run;
        stack_len++;
    }
    //perform actual band info encoding
    start = 0;
    for (i = stack_len - 1; i >= 0; i--) {
        put_bits(&s->pb, 4, stackcb[i]);
        count = stackrun[i];
        memset(sce->zeroes + win*16 + start, !stackcb[i], count);
        //XXX: memset when band_type is also uint8_t
        for (j = 0; j < count; j++) {
            sce->band_type[win*16 + start] =  stackcb[i];
            start++;
        }
        while (count >= run_esc) {
            put_bits(&s->pb, run_bits, run_esc);
            count -= run_esc;
        }
        put_bits(&s->pb, run_bits, count);
    }
}

static void codebook_trellis_rate(AACEncContext *s, SingleChannelElement *sce,
                                  int win, int group_len, const float lambda)
{
    BandCodingPath path[120][12];
    int w, swb, cb, start, size;
    int i, j;
    const int max_sfb  = sce->ics.max_sfb;
    const int run_bits = sce->ics.num_windows == 1 ? 5 : 3;
    const int run_esc  = (1 << run_bits) - 1;
    int idx, ppos, count;
    int stackrun[120], stackcb[120], stack_len;
    float next_minbits = INFINITY;
    int next_mincb = 0;

    start = win*128;
    for (cb = 0; cb < 12; cb++) {
        path[0][cb].cost     = run_bits+4;
        path[0][cb].prev_idx = -1;
        path[0][cb].run      = 0;
    }
    for (swb = 0; swb < max_sfb; swb++) {
        size = sce->ics.swb_sizes[swb];
        if (sce->zeroes[win*16 + swb]) {
            float cost_stay_here = path[swb][0].cost;
            float cost_get_here  = next_minbits + run_bits + 4;
            if (   run_value_bits[sce->ics.num_windows == 8][path[swb][0].run]
                != run_value_bits[sce->ics.num_windows == 8][path[swb][0].run+1])
                cost_stay_here += run_bits;
            if (cost_get_here < cost_stay_here) {
                path[swb+1][0].prev_idx = next_mincb;
                path[swb+1][0].cost     = cost_get_here;
                path[swb+1][0].run      = 1;
            } else {
                path[swb+1][0].prev_idx = 0;
                path[swb+1][0].cost     = cost_stay_here;
                path[swb+1][0].run      = path[swb][0].run + 1;
            }
            next_minbits = path[swb+1][0].cost;
            next_mincb = 0;
            for (cb = 1; cb < 12; cb++) {
                path[swb+1][cb].cost = 61450;
                path[swb+1][cb].prev_idx = -1;
                path[swb+1][cb].run = 0;
            }
        } else {
            float minbits = next_minbits;
            int mincb = next_mincb;
            int startcb = sce->band_type[win*16+swb];
            next_minbits = INFINITY;
            next_mincb = 0;
            for(w = 0; w < group_len; w++)
                ff_aac_quantize_band(s->qcoefs + w*36, sce->coeffs + start + w*128, sce->scoeffs + start + w*128,
                                     sce->sf_idx[(win+w)*16+swb], size);
            for (cb = 0; cb < startcb; cb++) {
                path[swb+1][cb].cost = 61450;
                path[swb+1][cb].prev_idx = -1;
                path[swb+1][cb].run = 0;
            }
            for (cb = startcb; cb < 12; cb++) {
                float cost_stay_here, cost_get_here;
                float bits = 0.0f;
                for (w = 0; w < group_len; w++) {
                    bits += sfb_size_hcb(s, s->qcoefs + w*36, size, cb);
                }
                cost_stay_here = path[swb][cb].cost + bits;
                cost_get_here  = minbits            + bits + run_bits + 4;
                if (   run_value_bits[sce->ics.num_windows == 8][path[swb][cb].run]
                    != run_value_bits[sce->ics.num_windows == 8][path[swb][cb].run+1])
                    cost_stay_here += run_bits;
                if (cost_get_here < cost_stay_here) {
                    path[swb+1][cb].prev_idx = mincb;
                    path[swb+1][cb].cost     = cost_get_here;
                    path[swb+1][cb].run      = 1;
                } else {
                    path[swb+1][cb].prev_idx = cb;
                    path[swb+1][cb].cost     = cost_stay_here;
                    path[swb+1][cb].run      = path[swb][cb].run + 1;
                }
                if (path[swb+1][cb].cost < next_minbits) {
                    next_minbits = path[swb+1][cb].cost;
                    next_mincb = cb;
                }
            }
        }
        start += sce->ics.swb_sizes[swb];
    }

    //convert resulting path from backward-linked list
    stack_len = 0;
    idx       = 0;
    for (cb = 1; cb < 12; cb++)
        if (path[max_sfb][cb].cost < path[max_sfb][idx].cost)
            idx = cb;
    ppos = max_sfb;
    while (ppos > 0) {
        assert(idx >= 0);
        cb = idx;
        stackrun[stack_len] = path[ppos][cb].run;
        stackcb [stack_len] = cb;
        idx = path[ppos-path[ppos][cb].run+1][cb].prev_idx;
        ppos -= path[ppos][cb].run;
        stack_len++;
    }
    //perform actual band info encoding
    start = 0;
    for (i = stack_len - 1; i >= 0; i--) {
        put_bits(&s->pb, 4, stackcb[i]);
        count = stackrun[i];
        memset(sce->zeroes + win*16 + start, !stackcb[i], count);
        //XXX: memset when band_type is also uint8_t
        for (j = 0; j < count; j++) {
            sce->band_type[win*16 + start] =  stackcb[i];
            start++;
        }
        while (count >= run_esc) {
            put_bits(&s->pb, run_bits, run_esc);
            count -= run_esc;
        }
        put_bits(&s->pb, run_bits, count);
    }
}

/** Return the minimum scalefactor where the quantized coef does not clip. */
static av_always_inline uint8_t coef2minsf(float coef) {
    return av_clip_uint8(log2f(coef)*4 - 69 + SCALE_ONE_POS - SCALE_DIV_512);
}

/** Return the maximum scalefactor where the quantized coef is not zero. */
static av_always_inline uint8_t coef2maxsf(float coef) {
    return av_clip_uint8(log2f(coef)*4 +  6 + SCALE_ONE_POS - SCALE_DIV_512);
}

typedef struct TrellisPath {
    float cost;
    int prev;
} TrellisPath;

#define TRELLIS_STAGES 121
#define TRELLIS_STATES (SCALE_MAX_DIFF+1)

static void search_for_quantizers_anmr(AVCodecContext *avctx, AACEncContext *s,
                                       SingleChannelElement *sce,
                                       const float lambda)
{
    int q, w, w2, g, start = 0;
    int i, j;
    int idx;
    TrellisPath paths[TRELLIS_STAGES][TRELLIS_STATES];
    int bandaddr[TRELLIS_STAGES];
    int minq;
    float mincost;
    float q0f = FLT_MAX, q1f = 0.0f, qnrgf = 0.0f;
    int q0, q1, qcnt = 0;

    for (i = 0; i < 1024; i++) {
        float t = fabsf(sce->coeffs[i]);
        if (t > 0.0f) {
            q0f = FFMIN(q0f, t);
            q1f = FFMAX(q1f, t);
            qnrgf += t*t;
            qcnt++;
        }
    }

    if (!qcnt) {
        memset(sce->sf_idx, 0, sizeof(sce->sf_idx));
        memset(sce->zeroes, 1, sizeof(sce->zeroes));
        return;
    }

    //minimum scalefactor index is when minimum nonzero coefficient after quantizing is not clipped
    q0 = coef2minsf(q0f);
    //maximum scalefactor index is when maximum coefficient after quantizing is still not zero
    q1 = coef2maxsf(q1f);
    //av_log(NULL, AV_LOG_ERROR, "q0 %d, q1 %d\n", q0, q1);
    if (q1 - q0 > 60) {
        int q0low  = q0;
        int q1high = q1;
        //minimum scalefactor index is when maximum nonzero coefficient after quantizing is not clipped
        int qnrg = av_clip_uint8(log2f(sqrtf(qnrgf/qcnt))*4 - 31 + SCALE_ONE_POS - SCALE_DIV_512);
        q1 = qnrg + 30;
        q0 = qnrg - 30;
        //av_log(NULL, AV_LOG_ERROR, "q0 %d, q1 %d\n", q0, q1);
        if (q0 < q0low) {
            q1 += q0low - q0;
            q0  = q0low;
        } else if (q1 > q1high) {
            q0 -= q1 - q1high;
            q1  = q1high;
        }
    }
    //av_log(NULL, AV_LOG_ERROR, "q0 %d, q1 %d\n", q0, q1);

    for (i = 0; i < TRELLIS_STATES; i++) {
        paths[0][i].cost    = 0.0f;
        paths[0][i].prev    = -1;
    }
    for (j = 1; j < TRELLIS_STAGES; j++) {
        for (i = 0; i < TRELLIS_STATES; i++) {
            paths[j][i].cost    = INFINITY;
            paths[j][i].prev    = -2;
        }
    }
    idx = 1;
    for (w = 0; w < sce->ics.num_windows; w += sce->ics.group_len[w]) {
        start = w*128;
        for (g = 0; g < sce->ics.num_swb; g++) {
            const float *coefs = sce->coeffs + start;
            float qmin, qmax;
            int nz = 0;

            bandaddr[idx] = w * 16 + g;
            qmin = INT_MAX;
            qmax = 0.0f;
            for (w2 = 0; w2 < sce->ics.group_len[w]; w2++) {
                FFPsyBand *band = &s->psy.ch[s->cur_channel].psy_bands[(w+w2)*16+g];
                if (band->energy <= band->threshold || band->threshold == 0.0f) {
                    sce->zeroes[(w+w2)*16+g] = 1;
                    continue;
                }
                sce->zeroes[(w+w2)*16+g] = 0;
                nz = 1;
                for (i = 0; i < sce->ics.swb_sizes[g]; i++) {
                    float t = fabsf(coefs[w2*128+i]);
                    if (t > 0.0f)
                        qmin = FFMIN(qmin, t);
                    qmax = FFMAX(qmax, t);
                }
            }
            if (nz) {
                int minscale, maxscale;
                float minrd = INFINITY;
                float maxval;
                //minimum scalefactor index is when minimum nonzero coefficient after quantizing is not clipped
                minscale = coef2minsf(qmin);
                //maximum scalefactor index is when maximum coefficient after quantizing is still not zero
                maxscale = coef2maxsf(qmax);
                minscale = av_clip(minscale - q0, 0, TRELLIS_STATES - 1);
                maxscale = av_clip(maxscale - q0, 0, TRELLIS_STATES);
                maxval = find_max_val(sce->ics.group_len[w], sce->ics.swb_sizes[g], sce->scoeffs + start);
                for (q = minscale; q < maxscale; q++) {
                    float dist = 0;
                    int cb = find_min_book(maxval, sce->sf_idx[w*16+g]);
                    for (w2 = 0; w2 < sce->ics.group_len[w]; w2++) {
                        FFPsyBand *band = &s->psy.ch[s->cur_channel].psy_bands[(w+w2)*16+g];
                        dist += sfb_rd_cost(s, coefs + w2*128, sce->scoeffs + start + w2*128, sce->ics.swb_sizes[g],
                                             q + q0, cb, lambda / band->threshold);
                    }
                    minrd = FFMIN(minrd, dist);

                    for (i = 0; i < q1 - q0; i++) {
                        float cost;
                        cost = paths[idx - 1][i].cost + dist
                               + ff_aac_scalefactor_bits[q - i + SCALE_DIFF_ZERO];
                        if (cost < paths[idx][q].cost) {
                            paths[idx][q].cost    = cost;
                            paths[idx][q].prev    = i;
                        }
                    }
                }
            } else {
                for (q = 0; q < q1 - q0; q++) {
                    paths[idx][q].cost = paths[idx - 1][q].cost + 1;
                    paths[idx][q].prev = q;
                }
            }
            sce->zeroes[w*16+g] = !nz;
            start += sce->ics.swb_sizes[g];
            idx++;
        }
    }
    idx--;
    mincost = paths[idx][0].cost;
    minq    = 0;
    for (i = 1; i < TRELLIS_STATES; i++) {
        if (paths[idx][i].cost < mincost) {
            mincost = paths[idx][i].cost;
            minq = i;
        }
    }
    while (idx) {
        sce->sf_idx[bandaddr[idx]] = minq + q0;
        minq = paths[idx][minq].prev;
        idx--;
    }
    //set the same quantizers inside window groups
    for (w = 0; w < sce->ics.num_windows; w += sce->ics.group_len[w])
        for (g = 0;  g < sce->ics.num_swb; g++)
            for (w2 = 1; w2 < sce->ics.group_len[w]; w2++)
                sce->sf_idx[(w+w2)*16+g] = sce->sf_idx[w*16+g];
}

/**
 * two-loop quantizers search taken from ISO 13818-7 Appendix C
 */
static void search_for_quantizers_twoloop(AVCodecContext *avctx,
                                          AACEncContext *s,
                                          SingleChannelElement *sce,
                                          const float lambda)
{
    int start = 0, i, w, w2, g;
    int destbits = avctx->bit_rate * 1024.0 / avctx->sample_rate / avctx->channels;
    float dists[128], uplims[128];
    float maxvals[128];
    int fflag, minscaler;
    int its  = 0;
    int allz = 0;
    float minthr = INFINITY;

    //XXX: some heuristic to determine initial quantizers will reduce search time
    memset(dists, 0, sizeof(dists));
    //determine zero bands and upper limits
    for (w = 0; w < sce->ics.num_windows; w += sce->ics.group_len[w]) {
        for (g = 0;  g < sce->ics.num_swb; g++) {
            int nz = 0;
            float uplim = 0.0f;
            for (w2 = 0; w2 < sce->ics.group_len[w]; w2++) {
                FFPsyBand *band = &s->psy.ch[s->cur_channel].psy_bands[(w+w2)*16+g];
                uplim += band->threshold;
                if (band->energy <= band->threshold || band->threshold == 0.0f) {
                    sce->zeroes[(w+w2)*16+g] = 1;
                    continue;
                }
                nz = 1;
            }
            uplims[w*16+g] = uplim *512;
            sce->zeroes[w*16+g] = !nz;
            if (nz)
                minthr = FFMIN(minthr, uplim);
            allz |= nz;
        }
    }
    for (w = 0; w < sce->ics.num_windows; w += sce->ics.group_len[w]) {
        for (g = 0;  g < sce->ics.num_swb; g++) {
            if (sce->zeroes[w*16+g]) {
                sce->sf_idx[w*16+g] = SCALE_ONE_POS;
                continue;
            }
            sce->sf_idx[w*16+g] = SCALE_ONE_POS + FFMIN(log2f(uplims[w*16+g]/minthr)*4,59);
        }
    }

    if (!allz)
        return;

    for (w = 0; w < sce->ics.num_windows; w += sce->ics.group_len[w]) {
        start = w*128;
        for (g = 0;  g < sce->ics.num_swb; g++) {
            const float *scaled = sce->scoeffs + start;
            maxvals[w*16+g] = find_max_val(sce->ics.group_len[w], sce->ics.swb_sizes[g], scaled);
            start += sce->ics.swb_sizes[g];
        }
    }

    //perform two-loop search
    //outer loop - improve quality
    do {
        int tbits, qstep;
        minscaler = sce->sf_idx[0];
        //inner loop - quantize spectrum to fit into given number of bits
        qstep = its ? 1 : 32;
        do {
            int prev = -1;
            tbits = 0;
            fflag = 0;
            for (w = 0; w < sce->ics.num_windows; w += sce->ics.group_len[w]) {
                start = w*128;
                for (g = 0;  g < sce->ics.num_swb; g++) {
                    const float *coefs = sce->coeffs + start;
                    const float *scaled = sce->scoeffs + start;
                    int bits = 0;
                    int cb;
                    float dist = 0.0f;

                    if (sce->zeroes[w*16+g] || sce->sf_idx[w*16+g] >= 218) {
                        start += sce->ics.swb_sizes[g];
                        continue;
                    }
                    minscaler = FFMIN(minscaler, sce->sf_idx[w*16+g]);
                    cb = find_min_book(maxvals[w*16+g], sce->sf_idx[w*16+g]);
                    for (w2 = 0; w2 < sce->ics.group_len[w]; w2++) {
                        dist += sfb_rd_cost(s, coefs + w2*128, scaled + w2*128, sce->ics.swb_sizes[g],
                                            sce->sf_idx[w*16+g], cb, 1.0f);
                        bits += sfb_size_hcb(s, s->qcoefs, sce->ics.swb_sizes[g], cb);
                    }
                    dists[w*16+g] = dist - bits;
                    if (prev != -1) {
                        bits += ff_aac_scalefactor_bits[sce->sf_idx[w*16+g] - prev + SCALE_DIFF_ZERO];
                    }
                    tbits += bits;
                    start += sce->ics.swb_sizes[g];
                    prev = sce->sf_idx[w*16+g];
                }
            }
            if (tbits > destbits) {
                for (i = 0; i < 128; i++)
                    if (sce->sf_idx[i] < 218 - qstep)
                        sce->sf_idx[i] += qstep;
            } else {
                for (i = 0; i < 128; i++)
                    if (sce->sf_idx[i] > 60 - qstep)
                        sce->sf_idx[i] -= qstep;
            }
            qstep >>= 1;
            if (!qstep && tbits > destbits*1.02 && sce->sf_idx[0] < 217)
                qstep = 1;
        } while (qstep);

        fflag = 0;
        minscaler = av_clip(minscaler, 60, 255 - SCALE_MAX_DIFF);
        for (w = 0; w < sce->ics.num_windows; w += sce->ics.group_len[w]) {
            for (g = 0; g < sce->ics.num_swb; g++) {
                int prevsc = sce->sf_idx[w*16+g];
                if (dists[w*16+g] > uplims[w*16+g] && sce->sf_idx[w*16+g] > 60) {
                    if (find_min_book(maxvals[w*16+g], sce->sf_idx[w*16+g]-1))
                        sce->sf_idx[w*16+g]--;
                    else //Try to make sure there is some energy in every band
                        sce->sf_idx[w*16+g]-=2;
                }
                sce->sf_idx[w*16+g] = av_clip(sce->sf_idx[w*16+g], minscaler, minscaler + SCALE_MAX_DIFF);
                sce->sf_idx[w*16+g] = FFMIN(sce->sf_idx[w*16+g], 219);
                if (sce->sf_idx[w*16+g] != prevsc)
                    fflag = 1;
                sce->band_type[w*16+g] = find_min_book(maxvals[w*16+g], sce->sf_idx[w*16+g]);
            }
        }
        its++;
    } while (fflag && its < 10);
}

static void search_for_quantizers_faac(AVCodecContext *avctx, AACEncContext *s,
                                       SingleChannelElement *sce,
                                       const float lambda)
{
    int start = 0, i, w, w2, g;
    float uplim[128], maxq[128];
    int minq, maxsf;
    float distfact = ((sce->ics.num_windows > 1) ? 85.80 : 147.84) / lambda;
    int last = 0, lastband = 0, curband = 0;
    float avg_energy = 0.0;
    if (sce->ics.num_windows == 1) {
        start = 0;
        for (i = 0; i < 1024; i++) {
            if (i - start >= sce->ics.swb_sizes[curband]) {
                start += sce->ics.swb_sizes[curband];
                curband++;
            }
            if (sce->coeffs[i]) {
                avg_energy += sce->coeffs[i] * sce->coeffs[i];
                last = i;
                lastband = curband;
            }
        }
    } else {
        for (w = 0; w < 8; w++) {
            const float *coeffs = sce->coeffs + w*128;
            start = 0;
            for (i = 0; i < 128; i++) {
                if (i - start >= sce->ics.swb_sizes[curband]) {
                    start += sce->ics.swb_sizes[curband];
                    curband++;
                }
                if (coeffs[i]) {
                    avg_energy += coeffs[i] * coeffs[i];
                    last = FFMAX(last, i);
                    lastband = FFMAX(lastband, curband);
                }
            }
        }
    }
    last++;
    avg_energy /= last;
    if (avg_energy == 0.0f) {
        for (i = 0; i < FF_ARRAY_ELEMS(sce->sf_idx); i++)
            sce->sf_idx[i] = SCALE_ONE_POS;
        return;
    }
    for (w = 0; w < sce->ics.num_windows; w += sce->ics.group_len[w]) {
        start = w*128;
        for (g = 0; g < sce->ics.num_swb; g++) {
            float *coefs   = sce->coeffs + start;
            const int size = sce->ics.swb_sizes[g];
            int start2 = start, end2 = start + size, peakpos = start;
            float maxval = -1, thr = 0.0f, t;
            maxq[w*16+g] = 0.0f;
            if (g > lastband) {
                maxq[w*16+g] = 0.0f;
                start += size;
                for (w2 = 0; w2 < sce->ics.group_len[w]; w2++)
                    memset(coefs + w2*128, 0, sizeof(coefs[0])*size);
                continue;
            }
            for (w2 = 0; w2 < sce->ics.group_len[w]; w2++) {
                for (i = 0; i < size; i++) {
                    float t = coefs[w2*128+i]*coefs[w2*128+i];
                    maxq[w*16+g] = FFMAX(maxq[w*16+g], fabsf(coefs[w2*128 + i]));
                    thr += t;
                    if (sce->ics.num_windows == 1 && maxval < t) {
                        maxval  = t;
                        peakpos = start+i;
                    }
                }
            }
            if (sce->ics.num_windows == 1) {
                start2 = FFMAX(peakpos - 2, start2);
                end2   = FFMIN(peakpos + 3, end2);
            } else {
                start2 -= start;
                end2   -= start;
            }
            start += size;
            thr = pow(thr / (avg_energy * (end2 - start2)), 0.3 + 0.1*(lastband - g) / lastband);
            t   = 1.0 - (1.0 * start2 / last);
            uplim[w*16+g] = distfact / (1.4 * thr + t*t*t + 0.075);
        }
    }
    memset(sce->sf_idx, 0, sizeof(sce->sf_idx));
    for (w = 0; w < sce->ics.num_windows; w += sce->ics.group_len[w]) {
        start = w*128;
        for (g = 0;  g < sce->ics.num_swb; g++) {
            const float *coefs  = sce->coeffs  + start;
            const float *scaled = sce->scoeffs + start;
            const int size      = sce->ics.swb_sizes[g];
            int scf, prev_scf, step;
            int min_scf = -1, max_scf = 256;
            float curdiff;
            if (maxq[w*16+g] < 21.544) {
                sce->zeroes[w*16+g] = 1;
                start += size;
                continue;
            }
            sce->zeroes[w*16+g] = 0;
            scf  = prev_scf = av_clip(SCALE_ONE_POS - SCALE_DIV_512 - log2f(1/maxq[w*16+g])*16/3, 60, 218);
            step = 16;
            for (;;) {
                float dist = 0.0f;
                int quant_max;

                for (w2 = 0; w2 < sce->ics.group_len[w]; w2++) {
                    dist += sfb_rd_cost(s, coefs + w2*128, scaled + w2*128, sce->ics.swb_sizes[g],
                                        scf, ESC_BT, lambda);
                    dist -= sfb_size_hcb(s, s->qcoefs, sce->ics.swb_sizes[g], ESC_BT);
                }
                dist *= 1.0f / 512.0f / lambda;
                quant_max = quant(maxq[w*16+g], ff_aac_pow2sf_tab[POW_SF2_ZERO - scf + SCALE_ONE_POS - SCALE_DIV_512]);
                if (quant_max >= 8191) { // too much, return to the previous quantizer
                    sce->sf_idx[w*16+g] = prev_scf;
                    break;
                }
                prev_scf = scf;
                curdiff = fabsf(dist - uplim[w*16+g]);
                if (curdiff <= 1.0f)
                    step = 0;
                else
                    step = log2f(curdiff);
                if (dist > uplim[w*16+g])
                    step = -step;
                scf += step;
                scf = av_clip_uint8(scf);
                step = scf - prev_scf;
                if (FFABS(step) <= 1 || (step > 0 && scf >= max_scf) || (step < 0 && scf <= min_scf)) {
                    sce->sf_idx[w*16+g] = av_clip(scf, min_scf, max_scf);
                    break;
                }
                if (step > 0)
                    min_scf = prev_scf;
                else
                    max_scf = prev_scf;
            }
            start += size;
        }
    }
    minq = sce->sf_idx[0] ? sce->sf_idx[0] : INT_MAX;
    for (i = 1; i < 128; i++) {
        if (!sce->sf_idx[i])
            sce->sf_idx[i] = sce->sf_idx[i-1];
        else
            minq = FFMIN(minq, sce->sf_idx[i]);
    }
    if (minq == INT_MAX)
        minq = 0;
    minq = FFMIN(minq, SCALE_MAX_POS);
    maxsf = FFMIN(minq + SCALE_MAX_DIFF, SCALE_MAX_POS);
    for (i = 126; i >= 0; i--) {
        if (!sce->sf_idx[i])
            sce->sf_idx[i] = sce->sf_idx[i+1];
        sce->sf_idx[i] = av_clip(sce->sf_idx[i], minq, maxsf);
    }
}

static void search_for_quantizers_fast(AVCodecContext *avctx, AACEncContext *s,
                                       SingleChannelElement *sce,
                                       const float lambda)
{
    int i, w, w2, g;
    int minq = 255;

    memset(sce->sf_idx, 0, sizeof(sce->sf_idx));
    for (w = 0; w < sce->ics.num_windows; w += sce->ics.group_len[w]) {
        for (g = 0; g < sce->ics.num_swb; g++) {
            for (w2 = 0; w2 < sce->ics.group_len[w]; w2++) {
                FFPsyBand *band = &s->psy.ch[s->cur_channel].psy_bands[(w+w2)*16+g];
                if (band->energy <= band->threshold) {
                    sce->sf_idx[(w+w2)*16+g] = 218;
                    sce->zeroes[(w+w2)*16+g] = 1;
                } else {
                    sce->sf_idx[(w+w2)*16+g] = av_clip(SCALE_ONE_POS - SCALE_DIV_512 + log2f(band->threshold), 80, 218);
                    sce->zeroes[(w+w2)*16+g] = 0;
                }
                minq = FFMIN(minq, sce->sf_idx[(w+w2)*16+g]);
            }
        }
    }
    for (i = 0; i < 128; i++) {
        sce->sf_idx[i] = 140;
        //av_clip(sce->sf_idx[i], minq, minq + SCALE_MAX_DIFF - 1);
    }
    //set the same quantizers inside window groups
    for (w = 0; w < sce->ics.num_windows; w += sce->ics.group_len[w])
        for (g = 0;  g < sce->ics.num_swb; g++)
            for (w2 = 1; w2 < sce->ics.group_len[w]; w2++)
                sce->sf_idx[(w+w2)*16+g] = sce->sf_idx[w*16+g];
}

static void search_for_ms(AACEncContext *s, ChannelElement *cpe,
                          const float lambda)
{
    int start = 0, i, w, w2, g;
    float M[128], S[128];
    SingleChannelElement *sce0 = &cpe->ch[0];
    SingleChannelElement *sce1 = &cpe->ch[1];
    float *L34 = sce0->scoeffs, *R34 = sce1->scoeffs, *M34 = sce0->ret, *S34 = sce1->ret;
    if (!cpe->common_window)
        return;
    for (w = 0; w < sce0->ics.num_windows; w += sce0->ics.group_len[w]) {
        for (g = 0;  g < sce0->ics.num_swb; g++) {
            if (!cpe->ch[0].zeroes[w*16+g] && !cpe->ch[1].zeroes[w*16+g]) {
                float dist1 = 0.0f, dist2 = 0.0f;
                for (w2 = 0; w2 < sce0->ics.group_len[w]; w2++) {
                    FFPsyBand *band0 = &s->psy.ch[s->cur_channel+0].psy_bands[(w+w2)*16+g];
                    FFPsyBand *band1 = &s->psy.ch[s->cur_channel+1].psy_bands[(w+w2)*16+g];
                    float minthr = FFMIN(band0->threshold, band1->threshold);
                    float maxthr = FFMAX(band0->threshold, band1->threshold);
                    for (i = 0; i < sce0->ics.swb_sizes[g]; i++) {
                        M[i] = (sce0->coeffs[start+w2*128+i]
                              + sce1->coeffs[start+w2*128+i]) * 0.5;
                        S[i] =  M[i]
                              - sce1->coeffs[start+w2*128+i];
                    }
                    abs_pow34_v(M34, M, sce0->ics.swb_sizes[g]);
                    abs_pow34_v(S34, S, sce0->ics.swb_sizes[g]);
                    dist1 += sfb_rd_cost(s, sce0->coeffs + start + w2*128, L34 + start + w2*128,
                                          sce0->ics.swb_sizes[g], sce0->sf_idx[(w+w2)*16+g],
                                          sce0->band_type[(w+w2)*16+g], lambda / band0->threshold);
                    dist1 += sfb_rd_cost(s, sce1->coeffs + start + w2*128, R34 + start + w2*128,
                                          sce1->ics.swb_sizes[g], sce1->sf_idx[(w+w2)*16+g],
                                          sce1->band_type[(w+w2)*16+g], lambda / band1->threshold);
                    dist2 += sfb_rd_cost(s, M, M34, sce0->ics.swb_sizes[g], sce0->sf_idx[(w+w2)*16+g],
                                          sce0->band_type[(w+w2)*16+g], lambda / maxthr);
                    dist2 += sfb_rd_cost(s, S, S34, sce1->ics.swb_sizes[g], sce1->sf_idx[(w+w2)*16+g],
                                          sce1->band_type[(w+w2)*16+g], lambda / minthr);
                }
                cpe->ms_mask[w*16+g] = dist2 < dist1;
            }
            start += sce0->ics.swb_sizes[g];
        }
    }
}

AACCoefficientsEncoder ff_aac_coders[] = {
    {
        search_for_quantizers_faac,
        encode_window_bands_info,
        search_for_ms,
    },
    {
        search_for_quantizers_anmr,
        encode_window_bands_info,
        search_for_ms,
    },
    {
        search_for_quantizers_twoloop,
        codebook_trellis_rate,
        search_for_ms,
    },
    {
        search_for_quantizers_fast,
        encode_window_bands_info,
        search_for_ms,
    },
};
