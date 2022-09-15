import _Exampledwarf
import f90wrap.runtime
import logging
import numpy as np

class Dwarf(f90wrap.runtime.FortranModule):

    @staticmethod
    def do_dwarf_call(numomp, nproma, nlev, ngptot, ngptotg, nblocks):
        """
        _Exampledwarf.cloudsc_driver_print(NUMOMP, NPROMA, NLEV, NGPTOT, NGPTOTG, NBLOCKS) 

        Parameters
        ----------
        numomp : int
        nproma : int
        nlev   : int
        ngptot : int
        ngptotg: int
        nblocks: int

        """
        _Exampledwarf.f90wrap_cloudsc_driver_print(numomp=numomp, nproma=nproma, nlev=nlev, ngptot=ngptot, ngptotg=ngptotg, nblocks=nblocks)
        
        
    @staticmethod
    def do_dwarf_full_call(numomp, nproma, nlev, ngptot, nblocks, ngptotg,
                           ptsphy,
                           pt, pq, 
             #             tendency_cml, tendency_loc, 
                           buffer_cml, buffer_loc, 
                           pap,      paph, 
                           plu,      plude,    pmfu,     pmfd, 
                           pa,       pclv,     psupsat,
                           pcovptot, 
                           pfplsl,   pfplsn,   pfhpsl,   pfhpsn):
        """
        _Exampledwarf.cloudsc_driver(NUMOMP, NPROMA, NLEV, NGPTOT, NGPTOTG, NBLOCKS) 

        Parameters
        ----------
        numomp : int
        nproma : int
        nlev   : int
        ngptot : int
        ngptotg: int
        nblocks: int

        """
        _Exampledwarf.f90wrap_cloudsc_driver_test(numomp=numomp, nproma=nproma, nlev=nlev, ngptot=ngptot,  ngpblks=nblocks, ngptotg=ngptotg,
                                                  ptsphy=ptsphy,
                                                   pt=pt, pq=pq,
                                                #   tendency_cml=tendency_cml, tendency_loc=tendency_loc,
                                                   buffer_cml=buffer_cml, buffer_loc=buffer_loc,
                                                   pap=pap,      paph=paph,
                                                   plu=plu,      plude=plude,    pmfu=pmfu,     pmfd=pmfd,
                                                   pa=pa,       pclv=pclv,     psupsat=psupsat,
                                                   pcovptot=pcovptot,
                                                   pfplsl=pfplsl,   pfplsn=pfplsn,   pfhpsl=pfhpsl,   pfhpsn=pfhpsn)
dwarf = Dwarf()
numomp=1
nproma=100
nlev=100
ngptot=100
ngptotg=100
nblocks=1
ndim=5
ptsphy=1.
pt        = np.zeros((nproma,nlev  ,nblocks), order='F')
pq        = np.zeros((nproma,nlev  ,nblocks), order='F')
pap       = np.zeros((nproma,nlev  ,nblocks), order='F')
paph      = np.zeros((nproma,nlev+1,nblocks), order='F')
plu       = np.zeros((nproma,nlev  ,nblocks), order='F')
plude     = np.zeros((nproma,nlev  ,nblocks), order='F')
pmfu      = np.zeros((nproma,nlev  ,nblocks), order='F')
pmfd      = np.zeros((nproma,nlev  ,nblocks), order='F')
pa        = np.zeros((nproma,nlev  ,nblocks), order='F')
pclv      = np.zeros((nproma,nlev  ,nblocks), order='F')
psupsat   = np.zeros((nproma,nlev  ,nblocks), order='F')
pcovptot  = np.zeros((nproma,nlev  ,nblocks), order='F')
pfplsl    = np.zeros((nproma,nlev+1,nblocks), order='F')
pfplsn    = np.zeros((nproma,nlev+1,nblocks), order='F')
pfhpsl    = np.zeros((nproma,nlev+1,nblocks), order='F')
pfhpsn    = np.zeros((nproma,nlev+1,nblocks), order='F')
#loc_T     = np.zeros((nproma,nlev  ,nblocks), order='F')
#loc_Q     = np.zeros((nproma,nlev  ,nblocks), order='F')
#loc_CLD   = np.zeros((nproma,nlev  ,ndim, nblocks), order='F')
#CML_T     = np.zeros((nproma,nlev  ,nblocks), order='F')
#CML_Q     = np.zeros((nproma,nlev  ,nblocks), order='F')
#CML_CLD   = np.zeros((nproma,nlev  ,ndim,nblocks), order='F')
buffer_cml     = np.zeros((nproma,nlev,3+ndim,nblocks), order='F')
buffer_loc     = np.zeros((nproma,nlev,3+ndim,nblocks), order='F')

dwarf.do_dwarf_call(numomp, nproma, nlev, ngptot, ngptotg, nblocks)
dwarf.do_dwarf_full_call(numomp, nproma, nlev, ngptot, nblocks, ngptotg,
                         ptsphy,
                         pt,pq,
                         buffer_cml,buffer_loc,
                         pap, paph,
                         plu, plude, pmfu, pmfd,
                         pa,pclv,psupsat,
                         pcovptot,
                         pfplsl, pfplsn, pfhpsl, pfhpsn)
