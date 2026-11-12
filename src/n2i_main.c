/*
 * n2i_main.c - hp_tg_g_F_n2i / ofct_srec_g_F_n2i: Noise2Inverse on the
 *              q????.img / .tif data sets (rhp_c.c), the same front end as
 *              his_srec_g_F_n2i (n2i.h) in front of the unchanged bodies,
 *              linked with their main() renamed (-DN2I):
 *
 *   -DN2I_HP_TG      hp_tg_g_F_n2i     HiPic Dr RC RA0 rec K
 *                                      HiPic Dr L1 C1 L2 C2 RA0 rec K
 *                    180-degree scan, hp_tg_main; views a multiple of K
 *   -DN2I_OFCT_SREC  ofct_srec_g_F_n2i HiPic Rc Oy rangeList Dr RA0 rec K
 *                                      HiPic L1 C1 L2 C2 Oy rangeList Dr RA0 rec K
 *                    360-degree offset scan, ofct_srec_main; views a
 *                    multiple of 2K (the 3-argument memory estimate of
 *                    ofct_srec is not offered here)
 *
 * The sample views are split into K (2 or 4) interleaved subsets in angle
 * order and each is reconstructed into rec/0 .. rec/K-1 with RA0 shifted
 * by the angle of the subset's first view; every subset is a child
 * process of this program with RHIS_VIEWS=k/K (rhp_views.h, honoured by
 * rhp_c.c as by rhis_hp.c).  Both bodies are compiled with
 * -DONLY_CT_VIEWS, as in his_srec, so that the child reads the same views
 * (those before the first blank line of output.log) as the front end.
 */

#include "n2i.h"

#if defined(N2I_HP_TG)
extern int	hp_tg_main(int argc,char **argv);
#define N2I_BODY	hp_tg_main
static const N2iSpec	N2iBody={"hp_tg","HiPic Dr RC RA0 rec","HiPic Dr L1 C1 L2 C2 RA0 rec",7,10,0,NULL};
#elif defined(N2I_OFCT_SREC)
extern int	ofct_srec_main(int argc,char **argv);
#define N2I_BODY	ofct_srec_main
static const N2iSpec	N2iBody={"ofct_srec","HiPic Rc Oy rangeList Dr RA0 rec","HiPic L1 C1 L2 C2 Oy rangeList Dr RA0 rec",9,12,1,NULL};
#else
#error "build with -DN2I_HP_TG or -DN2I_OFCT_SREC"
#endif

int	main(int argc,char **argv)
{
	if (getenv("RHIS_N2I_CHILD")==NULL) return N2iParent(argc,argv,&N2iBody);
	return N2I_BODY(argc,argv);
}
