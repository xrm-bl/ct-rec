/* ============================================================================
 * cbp_pipe.h  -  1 スライス先行パイプライン用の薄い層 (hp_tg / ofct_srec)
 * ----------------------------------------------------------------------------
 * 目的: スライス z の GPU 計算(BeginCBP .. EndCBP)の間に、ホスト側で
 *       スライス z+1 の準備(投影コピー、黒投影検出、リング除去)を行う。
 *       そのために投影バッファを 2 面持ち、PrepareCBP がどちらを読むかを
 *       SelectCBPProjection() で切り替える。
 *
 * GPU ビルド (USE_GPU): cbp.cu が提供する AllocCBPProjection /
 *       SelectCBPProjection / BeginCBP / EndCBP をそのまま使う。
 *       BeginCBP は H2D 転送とカーネル投入で戻り、EndCBP が完了を待つ。
 *
 * CPU ビルド: cbp_thread_*.c には Begin/End がないので、ここで同期的に
 *       エミュレートする(選択バッファを P へコピーしてから CBP() を呼ぶ)。
 *       重ね合わせは起きないが、呼び出し側のループ構造は GPU 版と同一に
 *       なるので、CPU 版での結果比較がそのままループの検証になる。
 *
 * 使い方:
 *   CBP_PIPE_INIT(P,N,M);              InitCBP の直後に 1 回
 *   P2=CBP_PIPE_ALLOC();               2 面目(連続 M*N、GPU 版は pinned)
 *   CBP_PIPE_SELECT(Pb); CBP_PIPE_BEGIN(dr,r0,t0);   Pb を入力に計算開始
 *   F=CBP_PIPE_END();                  結果 (Float **, N*N)
 *
 * 契約: CBP_PIPE_END() が返す F は「次に CBP_PIPE_END() を呼ぶまで」有効。
 *   途中で CBP_PIPE_BEGIN() を呼んでも F の中身は変わらない(GPU 版は
 *   EndCBP がホスト側 f へ D2H した結果で、次の EndCBP まで書き換えない。
 *   CPU 版は CBP() の内部バッファを END 時に専用バッファへ写して返す)。
 *   呼び出し側は BEGIN(z+1) を投げた後に F(z) を読んでよい。
 * ==========================================================================*/
#ifndef CBP_PIPE_H
#define CBP_PIPE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cbp.h"

#ifdef USE_GPU

EXTERN Float	**AllocCBPProjection(void);
EXTERN void	SelectCBPProjection(Float **pp);

#define CBP_PIPE_INIT(P,N,M)		((void)0)
#define CBP_PIPE_ALLOC()		AllocCBPProjection()
#define CBP_PIPE_SELECT(pp)		SelectCBPProjection(pp)
#define CBP_PIPE_BEGIN(dr,r0,t0)	BeginCBP(dr,r0,t0)
#define CBP_PIPE_END()			EndCBP()

#else	/* ---- CPU backend: synchronous emulation ---- */

static Float	**cbp_pipe_P=NULL,	/* cbp 側の投影バッファ(InitCBP の戻り値) */
		**cbp_pipe_sel=NULL,	/* 次の Begin が読むバッファ */
		**cbp_pipe_F=NULL,	/* CBP() が返した内部バッファ */
		**cbp_pipe_Fout=NULL;	/* END が返す写し(次の END まで有効) */
static int	cbp_pipe_N=0,cbp_pipe_M=0;

static void	cbp_pipe_init(Float **P,int N,int M)
{
	int	y;

	cbp_pipe_P=cbp_pipe_sel=P; cbp_pipe_N=N; cbp_pipe_M=M;
	if ((cbp_pipe_Fout=(Float **)malloc(sizeof(Float *)*(size_t)N))==NULL ||
	    (cbp_pipe_Fout[0]=(Float *)malloc(sizeof(Float)*(size_t)N*(size_t)N))==NULL) {
	    fputs("cbp_pipe: memory allocation error for the CBP result copy.\n",stderr);
	    exit(1);
	}
	for (y=1; y<N; y++) cbp_pipe_Fout[y]=cbp_pipe_Fout[y-1]+N;
}

static Float	**cbp_pipe_alloc(void)
{
	Float	**pp;
	int	m;

	if ((pp=(Float **)malloc(sizeof(Float *)*(size_t)cbp_pipe_M))==NULL ||
	    (pp[0]=(Float *)malloc(sizeof(Float)*(size_t)cbp_pipe_M*
						 (size_t)cbp_pipe_N))==NULL)
	    return NULL;
	for (m=1; m<cbp_pipe_M; m++) pp[m]=pp[m-1]+cbp_pipe_N;
	return pp;
}

static void	cbp_pipe_select(Float **pp)
{
	cbp_pipe_sel=(pp!=NULL)?pp:cbp_pipe_P;
}

static void	cbp_pipe_begin(double dr,double r0,double t0)
{
	if (cbp_pipe_sel!=cbp_pipe_P) {
	    int	m;

	    #pragma omp parallel for
	    for (m=0; m<cbp_pipe_M; m++)
		memcpy(cbp_pipe_P[m],cbp_pipe_sel[m],
		       sizeof(Float)*(size_t)cbp_pipe_N);
	}
	cbp_pipe_F=CBP(dr,r0,t0);
}

static Float	**cbp_pipe_end(void)
{
	int	y;

	/* CBP() の内部バッファは次の CBP() で上書きされるので写して返す */
	#pragma omp parallel for
	for (y=0; y<cbp_pipe_N; y++)
	    memcpy(cbp_pipe_Fout[y],cbp_pipe_F[y],
		   sizeof(Float)*(size_t)cbp_pipe_N);
	return cbp_pipe_Fout;
}

#define CBP_PIPE_INIT(P,N,M)		cbp_pipe_init(P,N,M)
#define CBP_PIPE_ALLOC()		cbp_pipe_alloc()
#define CBP_PIPE_SELECT(pp)		cbp_pipe_select(pp)
#define CBP_PIPE_BEGIN(dr,r0,t0)	cbp_pipe_begin(dr,r0,t0)
#define CBP_PIPE_END()			cbp_pipe_end()

#endif	/* USE_GPU */

#endif	/* CBP_PIPE_H */
