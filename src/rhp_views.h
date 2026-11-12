/*
 * rhp_views.h - RHIS_VIEWS=k/K (Noise2Inverse): keep every K-th sample
 *               view, starting at the k-th in angle order; the I0 views
 *               and the dark are untouched.  Shared by the two rhp.h
 *               readers, rhis_hp.c (a.his) and rhp_c.c (q????.img / .tif),
 *               which call SelectViews() right after the qsort of hp->OL
 *               in InitReadHiPic (sample views first, in angle order).
 *
 * The caller (n2i.h: his_srec_g_F_n2i, hp_tg_g_F_n2i, ofct_srec_g_F_n2i)
 * shifts RA0 by the angle of view k, so the K subsets share one geometry.
 * Nt must be a multiple of K, otherwise the equal-step assumption of the
 * back-projection (dt = pi/M) would not match the true step K*dtheta of
 * every subset.
 */
#ifndef RHP_VIEWS_H
#define RHP_VIEWS_H

static void	SelectViews(HiPic *hp)
{
	char	*env;
	int	k,K,n=0,i;

	if ((env=getenv("RHIS_VIEWS"))==NULL || *env=='\0') return;
	if (sscanf(env,"%d/%d",&k,&K)!=2 || K<1 || k<0 || k>=K) {
	    (void)fprintf(stderr,"RHIS_VIEWS : expected k/K with 0 <= k < K.\n"); exit(1);
	}
	if (hp->Nt%K!=0) {
	    (void)fprintf(stderr,"RHIS_VIEWS=%s: %d sample views are not a multiple of K=%d.\n",env,hp->Nt,K);
	    exit(1);
	}
	for (i=0; i<hp->Nt; i++) if (i%K==k) hp->OL[n++]=hp->OL[i];
	memmove(hp->OL+n,hp->OL+hp->Nt,sizeof(OutputLog)*(size_t)(hp->Nq+1-hp->Nt));
	(void)fprintf(stderr,"RHIS_VIEWS=%d/%d: views %d, %d, ... (%d of %d sample views)\n",
		      k,K,hp->OL[0].q,(n>1)?hp->OL[1].q:hp->OL[0].q,n,hp->Nt);
	hp->Nt=n;
}

#endif	/* RHP_VIEWS_H */
