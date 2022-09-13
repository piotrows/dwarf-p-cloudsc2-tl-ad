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
    def do_dwarf_full_call(numomp, nproma, nlev, ngptot, ngptotg, nblocks,
                           ptsphy,
                           pt, pq, 
                           tendency_cml, tendency_loc, 
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
        _Exampledwarf.f90wrap_cloudsc_driver_test(numomp=numomp, nproma=nproma, nlev=nlev, ngptot=ngptot, ngptotg=ngptotg, nblocks=nblocks,
                                                  ptsphy=ptsphy,
                                                   pt=pt, pq=pq,
                                                   tendency_cml=tendency_cml, tendency_loc=tendency_loc,
                                                   pap=pap,      paph=paph,
                                                   plu=plu,      plude=plude,    pmfu=pmfu,     pmfd=pmfd,
                                                   pa=pa,       pclv=pclv,     psupsat=psupsat,
                                                   pcovptot=pcovptot,
                                                   pfplsl=pfplsl,   pfplsn=pfplsn,   pfhpsl=pfhpsl,   pfhpsn=pfhpsn)
dwarf = Dwarf()
numomp=1
nproma=1
nlev=1
ngptot=1
ngptotg=1
nblocks=1
ndim=1
ptsphy=1.
pt        = np.zeros((ngptot,nlev  ,ndim), order='F')
pq        = np.zeros((ngptot,nlev  ,ndim), order='F')
pap       = np.zeros((ngptot,nlev  ,ndim), order='F')
paph      = np.zeros((ngptot,nlev+1,ndim), order='F')
plu       = np.zeros((ngptot,nlev  ,ndim), order='F')
plude     = np.zeros((ngptot,nlev  ,ndim), order='F')
pmfu      = np.zeros((ngptot,nlev  ,ndim), order='F')
pmfd      = np.zeros((ngptot,nlev  ,ndim), order='F')
pa        = np.zeros((ngptot,nlev  ,ndim), order='F')
pclv      = np.zeros((ngptot,nlev  ,ndim), order='F')
psupsat   = np.zeros((ngptot,nlev  ,ndim), order='F')
pcovptot  = np.zeros((ngptot,nlev  ,ndim), order='F')
pfplsl    = np.zeros((ngptot,nlev+1,ndim), order='F')
pfplsn    = np.zeros((ngptot,nlev+1,ndim), order='F')
pfhpsl    = np.zeros((ngptot,nlev+1,ndim), order='F')
pfhpsn    = np.zeros((ngptot,nlev+1,ndim), order='F')

dwarf.do_dwarf_call(numomp, nproma, nlev, ngptot, ngptotg, nblocks)
dwarf.do_dwarf_full_call(numomp, nproma, nlev, ngptot, ngptotg, nblocks,
                         ptsphy,
                         pt,pq,
                         pap, paph,
                         plu, plude, pmfu, pmfd,
                         pa,pclv,psupsat,
                         pcovptot,
                         pfplsl, pfplsn, pfhpsl, pfhpsn)
