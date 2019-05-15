/**
 @file
 rezer~ - v.008 - sampling/looping object by raja rez
 
 @ingroup    MSP
 */

#include "ext.h"
#include "ext_obex.h"
#include "z_dsp.h"
#include "ext_buffer.h"

#ifdef MAC_VERSION
#define FRCNLN inline __attribute__((always_inline))
#endif
#ifdef WIN_VERSION
#define FRCNLN __forceinline
#endif

typedef struct _rezer {
    t_pxobject obj;
    t_buffer_ref *bf_ref;
    t_symbol *bufname;
    t_double sr;
    t_ptr_int nchan;
    t_ptr_int bframes;
    t_ptr_int bnc;
    t_ptr_int fad;
    t_ptr_int rfad;
    t_ptr_int xfad;
    t_ptr_int rxfad;
    t_ptr_int nend;
    t_ptr_int nstart;
    t_ptr_int nrend;
    t_ptr_int nrstart;
    t_ptr_int playC;
    t_ptr_int recC;
    t_ptr_int recdir;
    t_ptr_int mode;
    t_bool wrap;
    t_bool rwrap;
    t_bool change;
    t_bool dirty;
    t_bool bufmod;
    t_bool intrupt;
    t_double phprev;
    t_double rfprev;
    t_double speed;
    t_double fade;
    t_double pos;
    t_double xpos;
    t_double rpos;
    t_double rxpos;
    t_double vend;
    t_double vstart;
    t_double start;
    t_double end;
    t_double dur;
    t_double rdur;
    t_double srscale;
} t_rezer;


void rezer_sperf64(t_rezer *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long vecfrmz, long flags, void *userparam);
void rezer_mperf64(t_rezer *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long vecfz, long flags, void *userparam);
void rezer_dsp64(t_rezer *x, t_object *dsp64, short *count, double samprate, long mxvecsize, long flgs);
void rezer_modset(t_rezer *x);
void rezer_set(t_rezer *x, t_symbol *s);
void *rezer_new(t_symbol *s, long n);
void rezer_free(t_rezer *x);
void rezer_snoop(t_rezer *x, long n);
t_max_err rezer_notify(t_rezer *x, t_symbol *s, t_symbol *msg, void *sndr, void *data);
void rezer_mode(t_rezer *x, long n);
void rezer_fade(t_rezer *x, t_double f);
void rezer_vend(t_rezer *x, t_double nd);
void rezer_vstart(t_rezer *x, t_double strt);
void rezer_assist(t_rezer *x, void *b, long m, long a, char *s);
void rezer_dblclick(t_rezer *x);
static t_class *rezer_class; static t_symbol *ps_nothing, *ps_buffer_modified;
                                                              //easing function(fade/crossfade)
static FRCNLN double eas_func_up(double y1, double ramp, long fad)
{ return y1*(0.5*(1.0-cos((1.0-(((double)fad)/ramp))*PI))); }

static FRCNLN double eas_func_dwn(double y1, double ramp, long fad)
{ return y1*(0.5*(1.0-cos((((double)fad)/ramp)*PI))); }

static FRCNLN double eas_rec_up(double y1, double rmp, double odb, long fad)
{ return y1*(1.0-((odb*0.5)*(1.0-cos((1.0-((double)fad/rmp))*PI)))); }

static FRCNLN double eas_rec_dwn(double y1, double rmp, double odb, long fad)
{ return y1*(1.0-((odb*0.5)*(1.0-cos(((double)fad/rmp)*PI)))); }

static FRCNLN void interp_index
(t_ptr_int p, t_ptr_int *ndxz,t_ptr_int *ndxb,t_ptr_int *ndxc, t_double spd,t_ptr_int frms)
{
    t_ptr_int dr = (spd>=0)?1:-1;
    *ndxz = p - dr; frms -= 1;                              //index-calc for cubic interpolation
    if(*ndxz < 0) *ndxz = frms + *ndxz; else if(*ndxz > frms) *ndxz = *ndxz - frms;
    
    *ndxb = p + dr;
    if(*ndxb < 0) *ndxb = frms + *ndxb; else if(*ndxb > frms) *ndxb = *ndxb - frms;
    
    *ndxc = *ndxb + dr;
    if(*ndxc < 0) *ndxc = frms + *ndxc; else if(*ndxc > frms) *ndxc = *ndxc - frms;
    return;
}

static FRCNLN double interp(double f, double z, double a, double b, double c)
{ return (!f) ? a : ((((0.5*(c-z) + 1.5*(a-b))*f + (z - 2.5*a + b + b - 0.5*c))*f + (0.5*(b-z)))*f + a); }

static FRCNLN void regzger(t_ptr_int frm, long vecfz, t_double f, t_double dur, t_double spd, t_double fd,
                           t_double vstr, t_double vnd, t_ptr_int *nstr, t_ptr_int *nnd,
                           t_double *str, t_double *nd, t_bool *wr)
{              //registration and management of changes to playback boundaries
    //....*str = 'virtual' start of the buffer(in frames); *nd = 'virtual' end of buffer(in frames);
    //....frm = total duration in 'frames' of the buffer~; vecfz = signal vector size('vector frames');
    //....f = starting phase input(0-1); dur = duration input(0-1);
    //....vstr = 'vstart' message(0-1 in terms of a virtual start point in buffer~'s overall length);
    //....vnd = 'vend' message(0-1 like 'vstart' but for a virtual end point);
    //....*nstr(nstart) = internal start time in frames, dependent on 4 things:
    //....................'virtual' start/end points('vstart' and 'vend' messages),
    //..............................and the incoming start 'phase' and 'duration' signal inputs
    //....*nnd(nend) = internal end time in frames, like *nstrt(nstart) mentioned above;
    //....*wr = wraparound flag(for when start/end makes play wrap around virtual bounds within buffer~)
    t_double vdur; f = (f<0) ? (f*-1) : ((f>=0.999998)?0.999998:f); frm -= 1;
    if(vstr<vnd){ *str=vstr*frm; *nd=vnd*frm; }
    else if(vstr>vnd){ *str=vnd*frm; *nd=vstr*frm; }else{ *str=0.; *nd=frm; }; vdur=*nd-*str;
    if(dur>1.)dur=1.;else if(dur<-1.)dur=-1.;
    dur=dur*vdur;    if(dur==0){ dur=vecfz; *nstr=(f*vdur)+*str; *nnd=*nstr+dur; }
    else if(dur<0){ *nnd=(f*vdur)+*str; *nstr=*nnd+dur; if(*nstr<0)*nstr=*nd+*nstr; }
    else{ *nstr=(f*vdur)+*str; *nnd=*nstr+dur; }
    if(*nnd>*nd)*nnd=*str+(*nnd-*nd-1); if(*nnd<*str)*nnd=*nd-(*str-*nnd-1);
    if(*nnd==*nstr) { *nstr=0.; *nnd=*nd; }else if(*nnd<*nstr) *wr=1; else *wr=0; return;
}

static FRCNLN void recgcer(t_ptr_int frm, long vecfz, t_double f, t_double dur,
                           t_ptr_int *nstr, t_ptr_int *nnd, t_bool *wr)
{              //registration and management of changes to recording boundaries
    //....only difference from above is no need for 'start' since that is always 0. when recording
    f = (f<0) ? (f*-1) : ((f>=0.999998)?0.999998:f); frm -= 1;
    if(dur>1.)dur=1.;else if(dur<0.)dur=0.; dur=dur*frm;
    if(dur==0){ dur=vecfz; *nstr=(f*frm); *nnd=*nstr+dur; }
    else if(dur<0)
    { *nnd=(f*frm); *nstr=*nnd+dur; if(*nstr<0)*nstr=frm+*nstr; }else{ *nstr=(f*frm); *nnd=*nstr+dur; }
    if(*nnd>frm)*nnd=*nnd-frm-1; if(*nnd<0)*nnd=frm+*nnd+1;
    if(*nnd==*nstr) { *nstr = 0.; *nnd = frm-1; }else if(*nnd<*nstr) *wr=1; else *wr=0; return;
}

static FRCNLN void recwaits(t_ptr_int rxfad, t_double recc, t_ptr_int frm, long vecfz, t_double f,
                                t_double dur, t_double pos, t_double fade, t_double spd,
                                t_double *rpos, t_ptr_int *rfad, t_ptr_int *nstr, t_ptr_int *nnd,
                                t_ptr_int *recC, t_ptr_int *rdir, t_bool *wr, t_bool *drty)
{
    if(recc!=0)
    {
        recgcer(frm,vecfz<<1,f,dur,nstr,nnd,wr);
        if(recc>0){ *rpos=*nstr;*rdir=1; }else{ *rpos=*nnd;*rdir=-1; } *drty=1;
    }                                                         *rfad=fade; *recC=recc; return;
}

void ext_main(void *r)
{
    t_class *c = class_new("rezer~", (method)rezer_new, (method)rezer_free, sizeof(t_rezer), 0L, A_SYM, A_DEFLONG, 0);
    
    class_addmethod(c, (method)rezer_dsp64, "dsp64", A_CANT, 0);
    class_addmethod(c, (method)rezer_set, "set", A_SYM, 0);
    class_addmethod(c, (method)rezer_mode, "mode", A_LONG, 0);
    class_addmethod(c, (method)rezer_fade, "fade", A_FLOAT, 0);
    class_addmethod(c, (method)rezer_vend, "vend", A_FLOAT, 0);
    class_addmethod(c, (method)rezer_vstart, "vstart", A_FLOAT, 0);
    class_addmethod(c, (method)rezer_assist, "assist", A_CANT, 0);
    class_addmethod(c, (method)rezer_snoop, "snoop", A_NOTHING, 0);
    class_addmethod(c, (method)rezer_dblclick, "dblclick", A_CANT, 0);
    class_addmethod(c, (method)rezer_notify, "notify", A_CANT, 0);
    class_dspinit(c); class_register(CLASS_BOX, c); rezer_class = c;
    ps_nothing = gensym(""); ps_buffer_modified = gensym("buffer_modified");
}

void rezer_sperf64(t_rezer *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long vecfz, long flags, void *userparam)
{   //........inlets: starting 'phase', end, playControl, recordingControl, speed, recording-Inputs-L/R
    //........outlets: left-output, right-output (always stereo) ..phase-output
    t_double    *phase = ins[0];
    t_double    *dr = ins[1];
    t_double    *plyCtl = ins[2];
    t_double    *speed = ins[3];
    t_double    *rphase = ins[4];
    t_double    *rdr = ins[5];
    t_double    *rcCtl = ins[6];
    t_double    *ovdb = ins[7];
    t_double    *rInL = ins[8];
    t_double    *rInR = ins[9];
    t_double    *outL = outs[0];
    t_double    *outR = outs[1];
    t_double    *outPh = outs[2];
    t_double    *outrPh = outs[3];
    t_buffer_obj *buffer = buffer_ref_getobject(x->bf_ref); t_float *tab;
    t_bool      bufmod, dirty, wrap, rwrap, v, change; t_ptr_int n = vecfz;
    t_ptr_int   mode, playC, recC, frames, bnc, fad, xfad, rfad, rxfad, rcdx, rxdx, rdir;
    t_ptr_int   indx, indxz, indxb, indxc, xndx, xndxz, xndxb, xndxc, nend, nstart, nrend, nrstart;
    t_double    f, rf, phprev, rfprev, dif, frac, fric, odb, playc, recc, fade, start, vstart;
    t_double    oL,oR,xL,xR,oP,rP, spd, dur, rdur, rcL, rcR, pos, xpos, rpos, rxpos, end, vend, scldsr;
    bufmod=x->bufmod; wrap=x->wrap; rwrap=x->rwrap; dirty=x->dirty; if(bufmod){rezer_modset(x);}
    tab=buffer_locksamples(buffer); if(!tab || x->obj.z_disabled)goto zero;
    dur=x->dur; rdur=x->rdur; phprev=x->phprev; rfprev=x->rfprev; fade=x->fade; spd=x->speed;
    vstart=x->vstart; vend=x->vend; start=x->start; end=x->end; nrstart=x->nrstart; nrend=x->nrend;
    pos=x->pos; xpos=x->xpos; rpos=x->rpos; rxpos=x->rxpos;  v=0;
    nstart=x->nstart; nend=x->nend; fad=x->fad; rfad=x->rfad; xfad=x->xfad; rxfad=x->rxfad;
    frames=x->bframes; bnc=x->bnc; mode=x->mode; playC=x->playC; recC=x->recC; rdir=x->recdir;
    scldsr=x->srscale; if(!mode)change=x->change;else change=1;
    
    //...........................................Perform Loop
    
    while (n--)
    {
        f=*phase++; dur=*dr++; playc=*plyCtl++; spd=*speed++; spd *= scldsr;
        rf=*rphase++; rdur=*rdr++; recc=*rcCtl++; odb = *ovdb++; rcL=*rInL++; rcR=*rInR++;
        //....................................Starting/Stopping
        if (playc != playC)
        {
                if(xfad<0) //<-once on, fades take priority
                {
                    if (playc != 0)
                    {
                        if(mode==3)
                        { if(((recC==0)||(recC==-2))||(recC==2)){ playC = playc; fad = fade; } }
                        else{ playC = playc; fad = fade; }
                        regzger(frames,vecfz,f,dur,spd,fade,vstart,vend,&nstart,&nend,&start,&end,&wrap);
                        if(spd>=0) pos = nstart; else pos = nend;  change = 0;
                    }else{ playC = playc; fad = fade; }
                }//..Once playing: 'phprev!=phas' detects phase change,..'xfad' = crossfade, 'fad' = fadein/out
        }   //............'change' = vstart/end registration(waits for loop-point/phase-change)
        if(playC)
        {
            if (((phprev != f) && (fad<0))&&(xfad<0))
            {
                if(!mode)
                {
                    regzger(frames,vecfz,f,dur,spd,fade,vstart,vend,&nstart,&nend,&start,&end,&wrap);
                    xpos = pos; if(spd>=0) pos = nstart; else pos = nend; xfad=fade; change = 0;
                }
            }
            else if(!wrap)
            {
                if (pos>nend)
                {
                    xfad=fade; xpos=pos;
                    if((mode)||(change))
                    {
                        regzger(frames,vecfz,f,dur,spd,fade,vstart,vend,&nstart,&nend,&start,&end,&wrap);
                        pos = nstart; xfad=fade;
                        if(mode==1)
                        {
                            if(recc!=recC)
                            { recwaits(rxfad, recc, frames, vecfz<<1, rf, rdur, pos, fade, spd, &rpos,
                                       &rfad, &nrstart, &nrend, &recC, &rdir, &rwrap, &dirty); }
                            else if(recC)
                            {
                                recgcer(frames,vecfz<<1,rf,rdur,&nrstart,&nrend,&rwrap);
                                rxpos=rpos; rpos=nrstart; rxfad=fade;
                            }
                        }
                    }else pos=nstart+(pos-nend-1); change=0;
                }
                else if (pos<nstart)
                {
                    xfad=fade; xpos=pos;
                    if((mode)||(change))
                    {
                        regzger(frames,vecfz,f,dur,spd,fade,vstart,vend,&nstart,&nend,&start,&end,&wrap);
                        pos = nend; xfad=fade;
                        if(mode==1)
                        {
                            if(recc!=recC)
                            { recwaits(rxfad, recc, frames, vecfz<<1, rf, rdur, pos, fade, spd, &rpos,
                                       &rfad, &nrstart, &nrend, &recC, &rdir, &rwrap, &dirty); }
                            else if(recC)
                            {
                                recgcer(frames,vecfz<<1,rf,rdur,&nrstart,&nrend,&rwrap);
                                rxpos=rpos; rpos = nrend; rxfad=fade;
                            }
                        }
                    }else pos=nend-(nstart-pos-1); change=0;
                }
            }
            else
            {
                if ((pos>nend)&&(pos<nstart))
                {
                    xpos = pos; xfad=fade;
                    if((mode)||(change))
                    {
                        regzger(frames,vecfz,f,dur,spd,fade,vstart,vend,&nstart,&nend,&start,&end,&wrap);
                        if(spd>=0){ pos=nstart; }else{ pos=nend; } xfad=fade;
                        if(mode==1)
                        {
                            if(recc!=recC)
                            { recwaits(rxfad, recc, frames, vecfz<<1, rf, rdur, pos, fade, spd, &rpos,
                                       &rfad, &nrstart, &nrend, &recC, &rdir, &rwrap, &dirty); }
                            else if(recC)
                            {
                                recgcer(frames,vecfz<<1,rf,rdur,&nrstart,&nrend,&rwrap);
                                rxpos=rpos; if(rdir>0){ rpos=nrstart; }else{ rpos=nrend; } rxfad=fade;
                            }
                        }
                    }else
                    { if(spd>=0){ pos=nstart+(pos-nend-1); }else{ pos=nend-(nstart-pos-1); } } change=0;
                }if(pos>end){ pos=start+(pos-end-1); }else if(pos<start){ pos=end-(start-pos-1); }
            }
            if(pos>=frames)pos-=frames; if(pos<0)pos+=frames;                           indx=trunc(pos);
            if(spd>0){ frac=pos-indx; }else if(spd<0){ frac=1.0-(pos-indx); }else frac=0.0; pos+=spd;
            interp_index(indx,&indxz,&indxb,&indxc,spd,frames);
            xL = interp(frac, tab[indxz*bnc], tab[indx*bnc], tab[indxb*bnc], tab[indxc*bnc]);
            xR = interp(frac, tab[indxz*bnc+1], tab[indx*bnc+1], tab[indxb*bnc+1], tab[indxc*bnc+1]);
        //............................Crossfades(happen at loop-points and phase changes)
            if(xfad>=0)
            {
                if(xpos>=frames) xpos -= frames; if(xpos<0) xpos += frames;             xndx=trunc(xpos);
                if(spd>0){fric=xpos-xndx;}else if(spd<0){fric=1.0-(xpos-xndx);}else fric=0.0; xpos+=spd;
                interp_index(xndx,&xndxz,&xndxb,&xndxc,spd,frames);
                oL = interp(fric, tab[xndxz*bnc], tab[xndx*bnc], tab[xndxb*bnc], tab[xndxc*bnc]);
                oR = interp(fric, tab[xndxz*bnc+1], tab[xndx*bnc+1], tab[xndxb*bnc+1], tab[xndxc*bnc+1]);
                oL = eas_func_up(xL, fade, xfad) + eas_func_dwn(oL, fade, xfad);
                oR = eas_func_up(xR, fade, xfad) + eas_func_dwn(oR, fade, xfad);
                xfad--;
            }                   //........................Fade-In
            else if(fad>=0){ oL=eas_func_up(xL,fade,fad); oR=eas_func_up(xR,fade,fad); fad--; }
            else{ oL=xL; oR=xR; } phprev=f;  //..<-Regular Playback + Phase-History
        }
        else
        {       //.......................................Fade-Out
            if(fad>=0)
            {
                if(pos>=frames) pos -= frames; if(pos<0) pos += frames;                 indx=trunc(pos);
                if(spd>0){frac=pos-indx;}else if(spd<0){frac=1.0-(pos-indx);}else frac=0.0; pos+=spd;
                interp_index(indx,&indxz,&indxb,&indxc,spd,frames);
                xL = interp(frac, tab[indxz*bnc], tab[indx*bnc], tab[indxb*bnc], tab[indxc*bnc]);
                xR = interp(frac, tab[indxz*bnc+1], tab[indx*bnc+1], tab[indxb*bnc+1], tab[indxc*bnc+1]);
                oL = eas_func_dwn(xL, fade, fad); oR = eas_func_dwn(xR, fade, fad);     fad--;
            }else{ pos = oL = oR = 0.; }  //........<-Everything Off
        }                                   oP = pos/(double)frames; //<-Sample-Index Converted To Phase
        //............................................................................................
        //..............................................RECORDING.....................................
        if((mode>=2)||(mode==0))
        {
            if(((recc!=recC) && (rfad<0))&&(rxfad<0))
            {
                if(mode==3)
                {
                    if(recc==0){ vend=rpos/(double)frames; }
                    else
                    {
                        if((recc==1)||(recc==-1)){ vstart = rf; }
                        if((recc==2)||(recc==-2))
                        {
                            vend=rpos/(double)frames;
                            recgcer(frames,vecfz<<1,vstart,vend-vstart,&nrstart,&nrend,&rwrap);
                            rxpos=rpos; rxfad=fade; if(recc>0){rpos=nrstart; rdir=1;}else{rpos=nrend; rdir=-1;}
                        }
                    } v=1;
                }
                if((recc<2)&&(recc>-2))
                {
                       recwaits(rxfad, recc, frames, vecfz<<1, rf, rdur, pos, fade, spd, &rpos,
                                 &rfad, &nrstart, &nrend, &recC, &rdir, &rwrap, &dirty);
                }
                else
                {
                    recwaits(rxfad, recc, frames, vecfz<<1, vstart, vend-vstart, pos, fade, spd, &rpos,
                             &rfad, &nrstart, &nrend, &recC, &rdir, &rwrap, &dirty);
                    regzger(frames,vecfz,f,dur,spd,fade,vstart,vend,&nstart,&nend,&start,&end,&wrap);
                }
            }
        }   //..Once recording: 'rfprev!=rf' detects phase change,..'rxfad'=crossfade, 'rfad'=fadein/out
        if(recC != 0)
        {
            if(((rfprev != rf) && (rfad<0))&&(rxfad<0))
            {
                if(!mode)
                {
                    recgcer(frames,vecfz<<1,rf,rdur,&nrstart,&nrend,&rwrap);
                    rxpos=rpos; rxfad=fade; if(recc>0){rpos=nrstart; rdir=1;}else{rpos=nrend; rdir=-1;}
                }
            }else if(!rwrap)
            {
                if(rpos>nrend)
                {
                        if(mode==2)recgcer(frames,vecfz<<1,rf,rdur,&nrstart,&nrend,&rwrap);
                        rxpos=rpos; rpos=nrstart; rxfad=fade;
                }else if(rpos<nrstart)
                {
                        if(mode==2)recgcer(frames,vecfz<<1,rf,rdur,&nrstart,&nrend,&rwrap);
                        rxpos=rpos; rpos=nrend; rxfad=fade;
                }
            }else
            {
                if((rpos>nrend)&&(rpos<nrstart))
                {
                    if(mode==2)recgcer(frames,vecfz<<1,rf,rdur,&nrstart,&nrend,&rwrap);
                        rxpos=rpos; rxfad=fade; if(rdir>0){ rpos=nrstart; }else{ rpos=nrend; }
                }
            }
            if(rpos>=frames) rpos-=frames; if(rpos<0) rpos+=frames;           rcdx=trunc(rpos); rpos+=rdir;
            //............................Crossfades(happen at loop-points and phase changes)
            if(rxfad>=0)
            {
                if(rxpos>=frames) rxpos-=frames; if(rxpos<0) rxpos+=frames;   rxdx=trunc(rxpos); rxpos+=rdir;
                tab[rxdx*bnc]=eas_func_dwn(rcL,fade,rxfad)+eas_rec_dwn(tab[rxdx*bnc],fade,1-odb,rxfad);
                tab[rxdx*bnc+1]=eas_func_dwn(rcR,fade,rxfad)+eas_rec_dwn(tab[rxdx*bnc+1],fade,1-odb,rxfad);
                tab[rcdx*bnc]=eas_func_up(rcL,fade,rxfad)+eas_rec_up(tab[rcdx*bnc],fade,1-odb,rxfad);
                tab[rcdx*bnc+1]=eas_func_up(rcR,fade,rxfad)+eas_rec_up(tab[rcdx*bnc+1],fade,1-odb,rxfad);
                rxfad--;
            }else if(rfad>=0)     //........................Fade-In
            {
                tab[rcdx*bnc]=eas_func_up(rcL,fade,rfad)+eas_rec_up(tab[rcdx*bnc],fade,1-odb,rfad);
                tab[rcdx*bnc+1]=eas_func_up(rcR,fade,rfad)+eas_rec_up(tab[rcdx*bnc+1],fade,1-odb,rfad); rfad--;
            }                                                     //....RegularRec + Rec-Phase-Histry....
            else{ tab[rcdx*bnc]=rcL+tab[rcdx*bnc]*odb; tab[rcdx*bnc+1]=rcR+tab[rcdx*bnc+1]*odb; };   rfprev=rf;
        }
        else
        {       //.......................................Fade-Out
            if(rfad>=0)
            {
                if(rpos>=frames) rpos-=frames; if(rpos<0) rpos+=frames;       rcdx=trunc(rpos); rpos+=rdir;
                tab[rcdx*bnc]=eas_func_dwn(rcL,fade,rfad)+eas_rec_dwn(tab[rcdx*bnc],fade,1-odb,rfad);
                tab[rcdx*bnc+1]=eas_func_dwn(rcR,fade,rfad)+eas_rec_dwn(tab[rcdx*bnc+1],fade,1-odb,rfad);
                rfad--; if(rfad==-1) dirty=0;
            }
        }                       rP = rpos/(double)frames; //.......<-Sample-Index Converted to Phase
        if(dirty)
        {
            if(spd!=rdir)
            { dif = fabs(rpos-pos);
                if(dif<=fade*2){ oL=eas_func_dwn(oL,fade*2,dif); oR=eas_func_dwn(oR,fade*2,dif); } }
        }
        *outL++ = oL; *outR++ = oR; *outPh++ = oP; *outrPh++ = rP;//.....<-Output
    }
    if(dirty)buffer_setdirty(buffer);           if(v){ x->vstart=vstart; x->vend=vend; }
    x->playC=playC; x->recC=recC; x->fad=fad; x->rfad=rfad; x->xfad=xfad; x->rxfad=rxfad;
    x->pos=pos; x->xpos=xpos; x->rpos=rpos; x->rxpos=rxpos; x->recdir=rdir; x->speed=spd;
    x->nstart=nstart; x->nend=nend; x->nrstart=nrstart; x->nrend=nrend; x->start=start; x->end=end;
    x->vstart=vstart; x->vend=vend;x->wrap=wrap; x->rwrap=rwrap; x->phprev=phprev; x->rfprev=rfprev;
    x->dur=dur; x->rdur=rdur; x->dirty=dirty; x->change=change; buffer_unlocksamples(buffer); return;
zero:
    while (n--) { *outL++ = 0.0; *outR++ = 0.0; *outPh++ = 0.0; *outrPh++ = 0.0; }
}

void rezer_modset(t_rezer *x)
{
    t_ptr_int bnc, frames; t_buffer_obj *b = buffer_ref_getobject(x->bf_ref);
    buffer_locksamples(b);
    if (b)
    {
        x->srscale=buffer_getsamplerate(b)/x->sr;
        frames=buffer_getframecount(b); bnc=buffer_getchannelcount(b);
        if((x->bframes != frames)||(x->bnc != bnc))
        {
            x->pos=x->xpos=x->rpos=x->rxpos=x->rfprev=x->phprev=0.;
            x->fad=x->rfad=x->xfad=x->rxfad=-1;  x->bframes=frames; x->bnc=bnc;
            x->nstart=0.; x->nend=frames-1; x->start=0.; x->end=frames-1; x->nrstart=0.; x->nrend=frames-1;
        }
    }
    buffer_unlocksamples(b); x->bufmod=0;
}
 
void rezer_set(t_rezer *x, t_symbol *s)
{
    x->bufname=s; if(!x->bf_ref)x->bf_ref=buffer_ref_new((t_object *)x, s);else buffer_ref_set(x->bf_ref, s);
    rezer_modset(x);
}

void rezer_mode(t_rezer *x, long n)
{ x->mode=n; post("rezer~ switched to mode %ld",n); if(n>0){ x->vstart=0.; x->vend=0.999999; x->change=1; } }

void rezer_fade(t_rezer *x, t_double f){ if(f<0)f=0-f; x->fade=f; }

void rezer_vend(t_rezer *x, t_double nd)
{ nd=(nd>1.) ? 1. : nd; nd=(nd<=0.) ? 0.002 : nd; if(x->mode!=3){x->vend=nd; x->change=1;} }

void rezer_vstart(t_rezer *x, t_double strt)
{ strt=(strt>=1.) ? 0.998 : strt; strt=(strt<=0.) ? 0. : strt; if(x->mode!=3){x->vstart=strt; x->change=1;} }

void rezer_dsp64(t_rezer *x, t_object *dsp64, short *count, double samprate, long mxvecsize, long flgs)
{
    x->sr = samprate;
    rezer_set(x, x->bufname);
    if(x->nchan > 1) { dsp_add64(dsp64, (t_object *)x, (t_perfroutine64)rezer_sperf64, 0, NULL); }
    else{ dsp_add64(dsp64, (t_object *)x, (t_perfroutine64)rezer_mperf64, 0, NULL); }
}

void rezer_dblclick(t_rezer *x) { buffer_view(buffer_ref_getobject(x->bf_ref)); }

void rezer_snoop(t_rezer *x, long n)
{
    post("recC: %ld", x->recC); post("playC: %ld", x->playC); post("pos: %.3f", x->pos);
    post("wrap: %ld", x->wrap); post("change: %ld", x->change);
    post("xfad: %ld", x->xfad); post("vstart: %.3f", x->vstart); post("vend: %.3f", x->vend);
    post("start: %.3f", x->start); post("end: %.3f", x->end); post("nend: %ld", x->nend);
    post("nstart: %ld", x->nstart); post("rpos: %.3f", x->rpos); post("rwrap: %ld", x->rwrap);
    post("rfad: %ld", x->rfad); post("rxfad: %ld", x->rxfad);
    post("nrend: %ld", x->nrend); post("nrstart: %ld", x->nrstart); post("intrupt: %ld", x->intrupt);
}

void rezer_assist(t_rezer *x, void *b, long m, long a, char *s)
{
    if (m == ASSIST_OUTLET) {switch (a)
        {
            case 0:    sprintf(s,"Left Out");    break;
            case 1:    if(x->nchan > 1) sprintf(s,"Right Out"); else sprintf(s,"Play Phase Out"); break;
            case 2:    if(x->nchan > 1) sprintf(s,"Play Phase Out"); else sprintf(s,"Rec Phase Out"); break;
            case 3:    sprintf(s, "Rec Phase Out"); break;
        }       }
    else {switch (a)
        {
            case 0:    sprintf(s,"(signal) Starting Phase/Index(0-1) (+ other messages)..."); break;
            case 1:    sprintf(s,"Duration Of Loop Window(0-1 as phase between virtual start/end"); break;
            case 2:    sprintf(s,"Play Control(0/1)");    break;
            case 3:    sprintf(s,"Speed(-∞/+∞)");    break;
            case 4:    sprintf(s,"Starting Recording-Phase/Index");    break;
            case 5:    sprintf(s,"Duration Of Rec Loop Window(0-1 as phase between real start/end"); break;
            case 6:    sprintf(s,"Recording Control(0/1)");    break;
            case 7:    sprintf(s,"Overdub Amount(Amplitude:0-1)");    break;
            case 8:    sprintf(s,"Recording Input Left");    break;
            case 9:    sprintf(s,"Recording Input Right");    break;
        }       }
}

void *rezer_new(t_symbol *s, long n)
{
    long chan = 2; if(n) chan = n;
    t_rezer *x = object_alloc(rezer_class);
    if(chan>1)
    {
        dsp_setup((t_pxobject *)x, 10);
        outlet_new((t_object *)x,"signal"); outlet_new((t_object *)x,"signal");
        outlet_new((t_object *)x,"signal"); outlet_new((t_object *)x,"signal");
    }
    else
    {
        dsp_setup((t_pxobject *)x, 9);
        outlet_new((t_object *)x,"signal"); outlet_new((t_object *)x,"signal");
        outlet_new((t_object *)x,"signal");
    }
    x->nchan = chan; x->bufname = s;
    x->bufmod = x->dirty = x->mode = x->wrap = x->rwrap = x->playC = x->recC = 0; x->fade = 512.;
    x->intrupt = x->change = x->phprev = x->rfprev = x->vstart = x->pos = x->xpos = x->rpos = x->rxpos = 0.;
    x->fad = x->xfad = x->rfad = x->rxfad= -1; x->recdir = x->vend = x->speed = 1.;
    post("rezer~ v.008 - mode %ld", x->mode); x->obj.z_misc |= Z_NO_INPLACE; return (x);
}

void rezer_free(t_rezer *x) { dsp_free((t_pxobject *)x); object_free(x->bf_ref); }

t_max_err rezer_notify(t_rezer *x, t_symbol *s, t_symbol *msg, void *sndr, void *data)
{ if(msg == ps_buffer_modified){ x->bufmod=1; } return buffer_ref_notify(x->bf_ref, s, msg, sndr, data); }

            /*_____________________________MONO PERF LOOP_______________________________________*/
void rezer_mperf64(t_rezer *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long vecfz, long flags, void *userparam)
{   //........inlets: starting 'phase', end, playControl, recordingControl, speed, recording-Inputs-L/R
    //........outlets: left-output, right-output (always stereo) ..phase-output
    t_double    *phase = ins[0];
    t_double    *dr = ins[1];
    t_double    *plyCtl = ins[2];
    t_double    *speed = ins[3];
    t_double    *rphase = ins[4];
    t_double    *rdr = ins[5];
    t_double    *rcCtl = ins[6];
    t_double    *ovdb = ins[7];
    t_double    *rInL = ins[8];
    t_double    *outL = outs[0];
    t_double    *outPh = outs[1];
    t_double    *outrPh = outs[2];
    t_buffer_obj *buffer = buffer_ref_getobject(x->bf_ref); t_float *tab;
    t_bool      bufmod, dirty, wrap, rwrap, v, change; t_ptr_int n = vecfz;
    t_ptr_int   mode, playC, recC, frames, bnc, fad, xfad, rfad, rxfad, rcdx, rxdx, rdir;
    t_ptr_int   indx, indxz, indxb, indxc, xndx, xndxz, xndxb, xndxc, nend, nstart, nrend, nrstart;
    t_double    f, rf, phprev, rfprev, dif, frac, fric, odb, playc, recc, fade, start, vstart;
    t_double    oL,xL,oP,rP, spd, dur, rdur, rcL, pos, xpos, rpos, rxpos, end, vend, scldsr;
    bufmod=x->bufmod; wrap=x->wrap; rwrap=x->rwrap; dirty=x->dirty; if(bufmod){rezer_modset(x);}
    tab=buffer_locksamples(buffer); if(!tab || x->obj.z_disabled)goto zero;
    dur=x->dur; rdur=x->rdur; phprev=x->phprev; rfprev=x->rfprev; fade=x->fade; spd=x->speed;
    vstart=x->vstart; vend=x->vend; start=x->start; end=x->end; nrstart=x->nrstart; nrend=x->nrend;
    pos=x->pos; xpos=x->xpos; rpos=x->rpos; rxpos=x->rxpos;  v=0;
    nstart=x->nstart; nend=x->nend; fad=x->fad; rfad=x->rfad; xfad=x->xfad; rxfad=x->rxfad;
    frames=x->bframes; bnc=x->bnc; mode=x->mode; playC=x->playC; recC=x->recC; rdir=x->recdir;
    scldsr=x->srscale; if(!mode)change=x->change;else change=1;
    
    //...........................................Perform Loop
    
    while (n--)
    {
        f=*phase++; dur=*dr++; playc=*plyCtl++; spd=*speed++; spd *= scldsr;
        rf=*rphase++; rdur=*rdr++; recc=*rcCtl++; odb = *ovdb++; rcL=*rInL++;
        //....................................Starting/Stopping
        if (playc != playC)
        {
            if(xfad<0) //<-once on, fades take priority
            {
                if (playc != 0)
                {
                    if(mode==3)
                    { if(((recC==0)||(recC==-2))||(recC==2)){ playC = playc; fad = fade; } }
                    else{ playC = playc; fad = fade; }
                    regzger(frames,vecfz,f,dur,spd,fade,vstart,vend,&nstart,&nend,&start,&end,&wrap);
                    if(spd>=0) pos = nstart; else pos = nend;  change = 0;
                }else{ playC = playc; fad = fade; }
            }//..Once playing: 'phprev!=phas' detects phase change,..'xfad' = crossfade, 'fad' = fadein/out
        }   //............'change' = vstart/end registration(waits for loop-point/phase-change)
        if(playC)
        {
            if (((phprev != f) && (fad<0))&&(xfad<0))
            {
                if(!mode)
                {
                    regzger(frames,vecfz,f,dur,spd,fade,vstart,vend,&nstart,&nend,&start,&end,&wrap);
                    xpos = pos; if(spd>=0) pos = nstart; else pos = nend; xfad=fade; change = 0;
                }
            }
            else if(!wrap)
            {
                if (pos>nend)
                {
                    xfad=fade; xpos=pos;
                    if((mode)||(change))
                    {
                        regzger(frames,vecfz,f,dur,spd,fade,vstart,vend,&nstart,&nend,&start,&end,&wrap);
                        pos = nstart; xfad=fade;
                        if(mode==1)
                        {
                            if(recc!=recC)
                            { recwaits(rxfad, recc, frames, vecfz<<1, rf, rdur, pos, fade, spd, &rpos,
                                       &rfad, &nrstart, &nrend, &recC, &rdir, &rwrap, &dirty); }
                            else if(recC)
                            {
                                recgcer(frames,vecfz<<1,rf,rdur,&nrstart,&nrend,&rwrap);
                                rxpos=rpos; rpos=nrstart; rxfad=fade;
                            }
                        }
                    }else pos=nstart+(pos-nend-1); change=0;
                }
                else if (pos<nstart)
                {
                    xfad=fade; xpos=pos;
                    if((mode)||(change))
                    {
                        regzger(frames,vecfz,f,dur,spd,fade,vstart,vend,&nstart,&nend,&start,&end,&wrap);
                        pos = nend; xfad=fade;
                        if(mode==1)
                        {
                            if(recc!=recC)
                            { recwaits(rxfad, recc, frames, vecfz<<1, rf, rdur, pos, fade, spd, &rpos,
                                       &rfad, &nrstart, &nrend, &recC, &rdir, &rwrap, &dirty); }
                            else if(recC)
                            {
                                recgcer(frames,vecfz<<1,rf,rdur,&nrstart,&nrend,&rwrap);
                                rxpos=rpos; rpos = nrend; rxfad=fade;
                            }
                        }
                    }else pos=nend-(nstart-pos-1); change=0;
                }
            }
            else
            {
                if ((pos>nend)&&(pos<nstart))
                {
                    xpos = pos; xfad=fade;
                    if((mode)||(change))
                    {
                        regzger(frames,vecfz,f,dur,spd,fade,vstart,vend,&nstart,&nend,&start,&end,&wrap);
                        if(spd>=0){ pos=nstart; }else{ pos=nend; } xfad=fade;
                        if(mode==1)
                        {
                            if(recc!=recC)
                            { recwaits(rxfad, recc, frames, vecfz<<1, rf, rdur, pos, fade, spd, &rpos,
                                       &rfad, &nrstart, &nrend, &recC, &rdir, &rwrap, &dirty); }
                            else if(recC)
                            {
                                recgcer(frames,vecfz<<1,rf,rdur,&nrstart,&nrend,&rwrap);
                                rxpos=rpos; if(rdir>0){ rpos=nrstart; }else{ rpos=nrend; } rxfad=fade;
                            }
                        }
                    }else
                    { if(spd>=0){ pos=nstart+(pos-nend-1); }else{ pos=nend-(nstart-pos-1); } } change=0;
                }if(pos>end){ pos=start+(pos-end-1); }else if(pos<start){ pos=end-(start-pos-1); }
            }
            if(pos>=frames)pos-=frames; if(pos<0)pos+=frames;                           indx=trunc(pos);
            if(spd>0){ frac=pos-indx; }else if(spd<0){ frac=1.0-(pos-indx); }else frac=0.0; pos+=spd;
            interp_index(indx,&indxz,&indxb,&indxc,spd,frames);
            xL = interp(frac, tab[indxz*bnc], tab[indx*bnc], tab[indxb*bnc], tab[indxc*bnc]);
            //............................Crossfades(happen at loop-points and phase changes)
            if(xfad>=0)
            {
                if(xpos>=frames) xpos -= frames; if(xpos<0) xpos += frames;             xndx=trunc(xpos);
                if(spd>0){fric=xpos-xndx;}else if(spd<0){fric=1.0-(xpos-xndx);}else fric=0.0; xpos+=spd;
                interp_index(xndx,&xndxz,&xndxb,&xndxc,spd,frames);
                oL = interp(fric, tab[xndxz*bnc], tab[xndx*bnc], tab[xndxb*bnc], tab[xndxc*bnc]);
                oL = eas_func_up(xL, fade, xfad) + eas_func_dwn(oL, fade, xfad);                xfad--;
            }                   //........................Fade-In
            else if(fad>=0){ oL=eas_func_up(xL,fade,fad); fad--; }else{ oL=xL; } phprev=f;
        }                                                             //..^Regular Playback + Phase-History
        else
        {       //.......................................Fade-Out
            if(fad>=0)
            {
                if(pos>=frames) pos -= frames; if(pos<0) pos += frames;                 indx=trunc(pos);
                if(spd>0){frac=pos-indx;}else if(spd<0){frac=1.0-(pos-indx);}else frac=0.0; pos+=spd;
                interp_index(indx,&indxz,&indxb,&indxc,spd,frames);
                xL = interp(frac, tab[indxz*bnc], tab[indx*bnc], tab[indxb*bnc], tab[indxc*bnc]);
                oL = eas_func_dwn(xL, fade, fad);                                                fad--;
            }else{ pos = oL = 0.; }  //........<-Everything Off
        }                                   oP = pos/(double)frames; //<-Sample-Index Converted To Phase
        //............................................................................................
        //..............................................RECORDING.....................................
        if((mode>=2)||(mode==0))
        {
            if(((recc!=recC) && (rfad<0))&&(rxfad<0))
            {
                if(mode==3)
                {
                    if(recc==0){ vend=rpos/(double)frames; }
                    else
                    {
                        if((recc==1)||(recc==-1)){ vstart = rf; }
                        if((recc==2)||(recc==-2))
                        {
                            vend=rpos/(double)frames;
                            recgcer(frames,vecfz<<1,vstart,vend-vstart,&nrstart,&nrend,&rwrap);
                            rxpos=rpos; rxfad=fade; if(recc>0){rpos=nrstart; rdir=1;}else{rpos=nrend; rdir=-1;}
                        }
                    } v=1;
                }
                if((recc<2)&&(recc>-2))
                {
                    recwaits(rxfad, recc, frames, vecfz<<1, rf, rdur, pos, fade, spd, &rpos,
                             &rfad, &nrstart, &nrend, &recC, &rdir, &rwrap, &dirty);
                }
                else
                {
                    recwaits(rxfad, recc, frames, vecfz<<1, vstart, vend-vstart, pos, fade, spd, &rpos,
                             &rfad, &nrstart, &nrend, &recC, &rdir, &rwrap, &dirty);
                    regzger(frames,vecfz,f,dur,spd,fade,vstart,vend,&nstart,&nend,&start,&end,&wrap);
                }
            }
        }   //..Once recording: 'rfprev!=rf' detects phase change,..'rxfad'=crossfade, 'rfad'=fadein/out
        if(recC != 0)
        {
            if(((rfprev != rf) && (rfad<0))&&(rxfad<0))
            {
                if(!mode)
                {
                    recgcer(frames,vecfz<<1,rf,rdur,&nrstart,&nrend,&rwrap);
                    rxpos=rpos; rxfad=fade; if(recc>0){rpos=nrstart; rdir=1;}else{rpos=nrend; rdir=-1;}
                }
            }else if(!rwrap)
            {
                if(rpos>nrend)
                {
                    if(mode==2)recgcer(frames,vecfz<<1,rf,rdur,&nrstart,&nrend,&rwrap);
                    rxpos=rpos; rpos=nrstart; rxfad=fade;
                }else if(rpos<nrstart)
                {
                    if(mode==2)recgcer(frames,vecfz<<1,rf,rdur,&nrstart,&nrend,&rwrap);
                    rxpos=rpos; rpos=nrend; rxfad=fade;
                }
            }else
            {
                if((rpos>nrend)&&(rpos<nrstart))
                {
                    if(mode==2)recgcer(frames,vecfz<<1,rf,rdur,&nrstart,&nrend,&rwrap);
                    rxpos=rpos; rxfad=fade; if(rdir>0){ rpos=nrstart; }else{ rpos=nrend; }
                }
            }
            if(rpos>=frames) rpos-=frames; if(rpos<0) rpos+=frames;           rcdx=trunc(rpos); rpos+=rdir;
            //............................Crossfades(happen at loop-points and phase changes)
            if(rxfad>=0)
            {
                if(rxpos>=frames) rxpos-=frames; if(rxpos<0) rxpos+=frames;   rxdx=trunc(rxpos); rxpos+=rdir;
                tab[rxdx*bnc]=eas_func_dwn(rcL,fade,rxfad)+eas_rec_dwn(tab[rxdx*bnc],fade,1-odb,rxfad);
                tab[rcdx*bnc]=eas_func_up(rcL,fade,rxfad)+eas_rec_up(tab[rcdx*bnc],fade,1-odb,rxfad); rxfad--;
            }else if(rfad>=0)     //........................Fade-In
            { tab[rcdx*bnc]=eas_func_up(rcL,fade,rfad)+eas_rec_up(tab[rcdx*bnc],fade,1-odb,rfad);   rfad--; }
            else{ tab[rcdx*bnc]=rcL+tab[rcdx*bnc]*odb; };   rfprev=rf; //..<-RegularRec + Rec-Phase-Histry
        }
        else
        {       //.......................................Fade-Out
            if(rfad>=0)
            {
                if(rpos>=frames) rpos-=frames; if(rpos<0) rpos+=frames;       rcdx=trunc(rpos); rpos+=rdir;
                tab[rcdx*bnc]=eas_func_dwn(rcL,fade,rfad)+eas_rec_dwn(tab[rcdx*bnc],fade,1-odb,rfad);
                rfad--; if(rfad==-1) dirty=0;
            }
        }                       rP = rpos/(double)frames; //.......<-Sample-Index Converted to Phase
        if(dirty)
        { if(spd!=rdir){ dif = fabs(rpos-pos); if(dif<=fade*2){ oL=eas_func_dwn(oL,fade*2,dif); } } }
        *outL++ = oL; *outPh++ = oP; *outrPh++ = rP;//.....<-Output
    }
    if(dirty)buffer_setdirty(buffer);           if(v){ x->vstart=vstart; x->vend=vend; }
    x->playC=playC; x->recC=recC; x->fad=fad; x->rfad=rfad; x->xfad=xfad; x->rxfad=rxfad;
    x->pos=pos; x->xpos=xpos; x->rpos=rpos; x->rxpos=rxpos; x->recdir=rdir; x->speed=spd;
    x->nstart=nstart; x->nend=nend; x->nrstart=nrstart; x->nrend=nrend; x->start=start; x->end=end;
    x->vstart=vstart; x->vend=vend;x->wrap=wrap; x->rwrap=rwrap; x->phprev=phprev; x->rfprev=rfprev;
    x->dur=dur; x->rdur=rdur; x->dirty=dirty; x->change=change; buffer_unlocksamples(buffer); return;
zero:
    while (n--) { *outL++ = 0.0; *outPh++ = 0.0; *outrPh++ = 0.0; }
}
