import _Exampledwarf
import f90wrap.runtime
import logging

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
        
        
dwarf = Dwarf()
numomp=1
nproma=1
nlev=1
ngptot=1
ngptotg=1
nblocks=1
dwarf.do_dwarf_call(numomp, nproma, nlev, ngptot, ngptotg, nblocks)
