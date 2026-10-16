/*
 * rhis.h - direct reader for a HiPic sequence file (a.his) plus the
 *          conv.bat / output.log that normally accompany it.
 *
 * The beamline writes every frame of a CT scan into one a.his file and
 * ships it together with conv.bat (how his2img / img_ave / ren / copy turn
 * the frames into dark.img and q????.img) and output.log (which q is an
 * I0 or a sample view, with time stamp and angle).  This module replays
 * conv.bat in memory, so that reconstruction programs can work on a.his
 * without ever creating the .img files:
 *
 *   his_rec  = ct_rec_c.c  built with -DHIS_INPUT, linked with rhis.c
 *   his_srec = hp_tg_ku.c  built with -DHIS_INPUT, linked with rhis.c
 *              (rhis.c also provides the rhp.h API: InitReadHiPic,
 *               ReadHiPic, ReadHiPicBand, TermReadHiPic)
 *
 * Only unsigned short frames (HIS type 2) are accepted.  Averages are
 * computed exactly like img_ave (integer sum, truncating division), so
 * dark / I0 data are bit-identical to the files conv.bat would produce.
 *
 * conv.bat grammar understood (case-insensitive command word):
 *   his2img <his>                      name of the his file (default a.his)
 *   img_ave <frame>... dark.img        dark = average of the listed frames
 *   img_ave <frame>... q<N>.img        q<N> = average of the listed frames
 *   ren|rename <frame> q<N>.img        q<N> = that frame
 *   copy <q<M>.img|frame> q<N>.img     q<N> = same data as q<M> (or frame)
 *   hp2do ... / rem / :: / @echo / echo / pause   ignored
 * <frame> is <letter><digits>.img as written by his2img (a0001.img is the
 * first frame of the file).  Any other line stops the program.
 *
 * Environment:
 *   RHIS_H   his file name (overrides the his2img line; default a.his)
 *   RHIS_C   conv.bat file name (default conv.bat)
 */
#ifndef RHIS_H
#define RHIS_H

#ifdef	__cplusplus
extern "C" {
#endif

typedef struct HisSet	HisSet;

/* open <dir>/a.his and interpret <dir>/conv.bat; exits on any error */
HisSet	*HisOpen(const char *dir);
void	HisClose(HisSet *s);

int	HisNx(const HisSet *s);		/* frame width  */
int	HisNy(const HisSet *s);		/* frame height */
int	HisNframes(const HisSet *s);	/* frames in the his file */
int	HisNq(const HisSet *s);		/* highest q number + 1 */
int	HisHasQ(const HisSet *s,int q);	/* q<N>.img is defined by conv.bat */
const char *HisPath(const HisSet *s);	/* path of the his file */

/* the image names on the hp2do line of conv.bat (D000 I000 E000 D180
   I180 E180), NULL where absent; returns how many were found (0 = no
   hp2do line).  Used by his_DO for the 180-degree axis estimate. */
int	HisHp2doNames(const HisSet *s,const char **names);
/* dark.img -> -1, q<N>.img defined by conv.bat -> N, otherwise -2 */
int	HisResolveName(const HisSet *s,const char *name);

/* rows y1..y2 (inclusive) of the dark frame / of q<q>, written to
   dst[(y-y1)*Nx + x].  Averaged entries are averaged like img_ave.
   Thread-safe for different q (shared read handle, positional reads). */
void	HisReadDarkBand(HisSet *s,int y1,int y2,unsigned short *dst);
void	HisReadQBand(HisSet *s,int q,int y1,int y2,unsigned short *dst);

#ifdef	__cplusplus
}
#endif

#endif	/* RHIS_H */
