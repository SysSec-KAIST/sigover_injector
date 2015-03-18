/**
 *
 * \section COPYRIGHT
 *
 * Copyright 2013-2014 The libLTE Developers. See the
 * COPYRIGHT file at the top-level directory of this distribution.
 *
 * \section LICENSE
 *
 * This file is part of the libLTE library.
 *
 * libLTE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * libLTE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * A copy of the GNU Lesser General Public License can be found in
 * the LICENSE file in the top-level directory of this distribution
 * and at http://www.gnu.org/licenses/.
 *
 */

#include <math.h>
#include <string.h>
#include "srslte/resampling/resample_arb.h"
#include "srslte/utils/debug.h"

float srslte_resample_arb_polyfilt[SRSLTE_RESAMPLE_ARB_N][SRSLTE_RESAMPLE_ARB_M] =
{{0,0.002400347599485495,-0.006922416132556366,0.0179104136912176,0.99453086623794,-0.008521087756729117,0.0008598969867484128,0.0004992625165376107},
{-0.001903604727400391,0.004479591950094871,-0.01525319260830623,0.04647449496926549,0.9910477342662829,-0.03275243420114668,0.008048813755373533,-0.001216900416836847},
{-0.001750442300940216,0.006728826416921727,-0.02407540632178267,0.07708575473589654,0.9841056525667189,-0.05473739187922162,0.01460652754040275,-0.002745266140572769},
{-0.001807302702047332,0.009130494591071001,-0.03332524241797466,0.1096377743266821,0.9737536341125557,-0.07444712408657775,0.0205085154280773,-0.004077596041064956},
{-0.002005048395707184,0.01166717081713966,-0.04292591596219218,0.1440068251440274,0.9600641782991848,-0.09187282739868929,0.02573649090563353,-0.005218582609802016},
{-0.002303606593778858,0.01431549924598744,-0.05279365431190471,0.1800496557331084,0.9431333663677619,-0.107022659158021,0.03028295984887362,-0.006178495199082676},
{-0.002673824939871474,0.01705308556587308,-0.06283262272105428,0.2176066284621515,0.9230793648715676,-0.1199244901666642,0.03414559686566755,-0.006951034565279722},
{-0.003099294850834212,0.01985022208215894,-0.07294100481214681,0.2564996326404178,0.9000412628278351,-0.1306207218667012,0.03733448898243594,-0.007554128492547957},
{-0.003564415127592977,0.02267702927238637,-0.08300660359189582,0.2965368855416621,0.874178327911914,-0.1391710078036285,0.03986356581130199,-0.007992692892036229},
{-0.004060222849803048,0.02549633298710617,-0.09291211362877161,0.3375102453232687,0.8456688916962831,-0.1456479420005528,0.04175429631970297,-0.00827816851803391},
{-0.004574770330592958,0.02827133859849648,-0.1025308044303699,0.379200979038518,0.8147072909915275,-0.1501397816831265,0.04303284689167361,-0.008423951215857236},
{-0.005100453570725489,0.0309593063669097,-0.1117325201826136,0.4213772298764557,0.7815036335229155,-0.1527444636485677,0.0437349357738538,-0.008441129184587313},
{-0.005625340687997995,0.03351877124912608,-0.1203802123857154,0.463798444519773,0.7462810500374923,-0.1535722712206338,0.04389261707032224,-0.008345136869130621},
{-0.006140888413286479,0.03590276234084221,-0.1283350405000867,0.5062161107254219,0.7092744343987509,-0.152741400055636,0.04355110899623683,-0.008147838953503964},
{-0.006634012711933725,0.03806645580467819,-0.1354549138065051,0.5483763739419955,0.6707272107491931,-0.1503798528838644,0.0427502160096194,-0.007865132034489236},
{-0.007094909626785393,0.03996043743797257,-0.1415969539549883,0.5900207719663111,0.6308907308946271,-0.1466186670140106,0.04153829696698895,-0.007508971112586246},
{-0.007508971112586246,0.04153829696698895,-0.1466186670140106,0.6308907308946271,0.5900207719663111,-0.1415969539549883,0.03996043743797257,-0.007094909626785393},
{-0.007865132034489236,0.0427502160096194,-0.1503798528838644,0.6707272107491931,0.5483763739419955,-0.1354549138065051,0.03806645580467819,-0.006634012711933725},
{-0.008147838953503964,0.04355110899623683,-0.152741400055636,0.7092744343987509,0.5062161107254219,-0.1283350405000867,0.03590276234084221,-0.006140888413286479},
{-0.008345136869130621,0.04389261707032224,-0.1535722712206338,0.7462810500374923,0.463798444519773,-0.1203802123857154,0.03351877124912608,-0.005625340687997995},
{-0.008441129184587313,0.0437349357738538,-0.1527444636485677,0.7815036335229155,0.4213772298764557,-0.1117325201826136,0.0309593063669097,-0.005100453570725489},
{-0.008423951215857236,0.04303284689167361,-0.1501397816831265,0.8147072909915275,0.379200979038518,-0.1025308044303699,0.02827133859849648,-0.004574770330592958},
{-0.00827816851803391,0.04175429631970297,-0.1456479420005528,0.8456688916962831,0.3375102453232687,-0.09291211362877161,0.02549633298710617,-0.004060222849803048},
{-0.007992692892036229,0.03986356581130199,-0.1391710078036285,0.874178327911914,0.2965368855416621,-0.08300660359189582,0.02267702927238637,-0.003564415127592977},
{-0.007554128492547957,0.03733448898243594,-0.1306207218667012,0.9000412628278351,0.2564996326404178,-0.07294100481214681,0.01985022208215894,-0.003099294850834212},
{-0.006951034565279722,0.03414559686566755,-0.1199244901666642,0.9230793648715676,0.2176066284621515,-0.06283262272105428,0.01705308556587308,-0.002673824939871474},
{-0.006178495199082676,0.03028295984887362,-0.107022659158021,0.9431333663677619,0.1800496557331084,-0.05279365431190471,0.01431549924598744,-0.002303606593778858},
{-0.005218582609802016,0.02573649090563353,-0.09187282739868929,0.9600641782991848,0.1440068251440274,-0.04292591596219218,0.01166717081713966,-0.002005048395707184},
{-0.004077596041064956,0.0205085154280773,-0.07444712408657775,0.9737536341125557,0.1096377743266821,-0.03332524241797466,0.009130494591071001,-0.001807302702047332},
{-0.002745266140572769,0.01460652754040275,-0.05473739187922162,0.9841056525667189,0.07708575473589654,-0.02407540632178267,0.006728826416921727,-0.001750442300940216},
{-0.001216900416836847,0.008048813755373533,-0.03275243420114668,0.9910477342662829,0.04647449496926549,-0.01525319260830623,0.004479591950094871,-0.001903604727400391},
{0.0004992625165376107,0.0008598969867484128,-0.008521087756729117,0.99453086623794,0.0179104136912176,-0.006922416132556366,0.002400347599485495,0}};


// TODO: use lte/utils/vector.h and Volk
cf_t srslte_resample_arb_dot_prod(cf_t* x, float *y, int len)
{
  cf_t res = 0+0*I;
  for(int i=0;i<len;i++){
    res += x[i]*y[i];
  }
  return res;
}

// Right-shift our window of samples
void srslte_resample_arb_push(srslte_resample_arb_t *q, cf_t x)
{
  memmove(&q->reg[1], &q->reg[0], (SRSLTE_RESAMPLE_ARB_M-1)*sizeof(cf_t));
  q->reg[0] = x;
}

// Initialize our struct
void srslte_resample_arb_init(srslte_resample_arb_t *q, float rate){
  memset(q->reg, 0, SRSLTE_RESAMPLE_ARB_M*sizeof(cf_t));
  q->acc = 0.0;
  q->rate = rate;
  q->step = (1/rate)*SRSLTE_RESAMPLE_ARB_N;
}

// Resample a block of input data
int srslte_resample_arb_compute(srslte_resample_arb_t *q, cf_t *input, cf_t *output, int n_in){
  int cnt = 0;
  int n_out = 0;
  int idx = 0;

  while(cnt < n_in)
  {
    *output = srslte_resample_arb_dot_prod(q->reg, srslte_resample_arb_polyfilt[idx], SRSLTE_RESAMPLE_ARB_M);
    output++;
    n_out++;
    q->acc += q->step;
    idx = (int)roundf(q->acc);
    while(idx >= SRSLTE_RESAMPLE_ARB_N){
      q->acc -= SRSLTE_RESAMPLE_ARB_N;
      idx -= SRSLTE_RESAMPLE_ARB_N;
      if(cnt < n_in)
        srslte_resample_arb_push(q, input[cnt++]);
    }
  }
  return n_out;
}