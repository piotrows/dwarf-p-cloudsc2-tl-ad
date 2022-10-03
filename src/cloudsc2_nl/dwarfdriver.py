import _Exampledwarf
import f90wrap.runtime
import logging
import numpy as np
from pathlib import Path
from collections import OrderedDict
from cloudsc2_inputs import load_input_parameters, load_input_fields, load_reference_fields
from cloudsc2_driver import arguments_from_fields 
class Dwarf(f90wrap.runtime.FortranModule):

    @staticmethod
    def do_dwarf_call_full(numomp, nproma, nlev, ngptot, ngptotg, nblocks,ptsphy):
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
        _Exampledwarf.f90wrap_cloudsc_driver_full_forpy(ptsphy=ptsphy, numomp=numomp, 
                                                        nproma=nproma, nlev=nlev,
                                                        ngptot=ngptot, ngpblks=nblocks,ngptotg=ngptotg)

    @staticmethod
    def do_dwarf_init_call( numomp, nproma, nlev, ngptot, nblocks, ngptotg,
                           ptsphy,
                           pt, pq, 
                           buffer_cml, buffer_loc, 
                           pap,      paph, 
                           plu,      plude,    pmfu,     pmfd, 
                           pa,       pclv,     psupsat,
                           pcovptot): 
        _Exampledwarf.f90wrap_cloudsc_driver_init( ptsphy=ptsphy,
                                                    numomp=numomp,
                                                    nproma=nproma, 
                                                  nlev=nlev,
                                                   ngptot=ngptot, 
                                                  ngpblks=nblocks, 
                                                  ngptotg=ngptotg,
                                                  output_pt=pt,
                                                  output_pq=pq,
                                                  output_buffer_cml=buffer_cml, 
                                                  output_buffer_loc=buffer_loc,
                                                  output_pap=pap,  
                                                  output_paph=paph,
                                                  output_plu=plu,
                                                  output_plude=plude,
                                                  output_pmfu=pmfu,
                                                  output_pmfd=pmfd,
                                                  output_pa=pa, 
                                                  output_pclv=pclv,
                                                  output_psupsat=psupsat,
                                                  output_pcovptot=pcovptot
                                                  )


    @staticmethod
    def do_dwarf_inittest_call( numomp, nproma, nlev, ngptot, nblocks, ngptotg,
                           ptsphy,
                           pt, pq, 
                           buffer_cml, buffer_loc, 
                           pap,      paph, 
                           plu,      plude,    pmfu,     pmfd, 
                           pa,       pclv,     psupsat,
                           pcovptot,
                           pfplsl,   pfplsn,   pfhpsl,   pfhpsn):

        _Exampledwarf.f90wrap_cloudsc_driver_inittest(
                           ptsphy=ptsphy,
                             numomp=numomp,
                             nproma=nproma, 
                           nlev=nlev,
                            ngptot=ngptot, 
                           ngpblks=nblocks, 
                           ngptotg=ngptotg,
                           output_pt=pt,
                           output_pq=pq,
                           output_buffer_cml=buffer_cml, 
                           output_buffer_loc=buffer_loc,
                           output_pap=pap,  
                           output_paph=paph,
                           output_plu=plu,
                           output_plude=plude,
                           output_pmfu=pmfu,
                           output_pmfd=pmfd,
                           output_pa=pa, 
                           output_pclv=pclv,
                           output_psupsat=psupsat,
                           output_pcovptot=pcovptot,
                           output_pfplsl=pfplsl,
                           output_pfplsn=pfplsn,
                           output_pfhpsl=pfhpsl,
                           output_pfhpsn=pfhpsn,
)


    @staticmethod
    def do_dwarf_validate_call(
                           numomp, nproma, nlev, ngptot, nblocks, ngptotg,
                           ptsphy,
                           pt, pq, 
                           buffer_cml, buffer_loc, 
                           pap,      paph, 
                           plu,      plude,    pmfu,     pmfd, 
                           pa,       pclv,     psupsat,
                           pcovptot, 
                           pfplsl,   pfplsn,   pfhpsl,   pfhpsn):
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
        _Exampledwarf.f90wrap_cloudsc_driver_validate( numomp=numomp, nproma=nproma, nlev=nlev, ngptot=ngptot, ngpblks=nblocks, ngptotg=ngptotg,
                                                 ptsphy=ptsphy,
                                                   input_pt=pt, input_pq=pq,
                                                   input_buffer_cml=buffer_cml, input_buffer_loc=buffer_loc,
                                                   input_pap=pap,      input_paph=paph,
                                                   input_plu=plu,      input_plude=plude,    input_pmfu=pmfu,     input_pmfd=pmfd,
                                                   input_pa=pa,        input_pclv=pclv,     input_psupsat=psupsat,
                                                   input_pcovptot=pcovptot, input_pfplsl=pfplsl,   input_pfplsn=pfplsn,   input_pfhpsl=pfhpsl,   input_pfhpsn=pfhpsn)

    @staticmethod
    def do_dwarf_full_call(
                           numomp, nproma, nlev, ngptot, nblocks, ngptotg,
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
        _Exampledwarf.f90wrap_cloudsc_driver_test(
                                                  numomp=numomp, nproma=nproma, nlev=nlev, ngptot=ngptot,  ngpblks=nblocks, ngptotg=ngptotg,
                                                  ptsphy=ptsphy,
                                                   pt=pt, pq=pq,
                                                #   tendency_cml=tendency_cml, tendency_loc=tendency_loc,
                                                   buffer_cml=buffer_cml, buffer_loc=buffer_loc,
                                                   pap=pap,      paph=paph,
                                                   plu=plu,      plude=plude,    pmfu=pmfu,     pmfd=pmfd,
                                                   pa=pa,       pclv=pclv,     psupsat=psupsat,
                                                   pcovptot=pcovptot,
                                                   pfplsl=pfplsl,   pfplsn=pfplsn,   pfhpsl=pfhpsl,   pfhpsn=pfhpsn)
    def examine_ndarray_flags(varname,ndvar):
        print ("Checking flags of array: ",varname)
        print (ndvar.flags)
        print ("End of flags of array: ",varname)
 

dwarf = Dwarf()
numomp=1
nproma=16384
nlev=137
ngptot=16384
ngptotg=16384
nblocks=1
ndim=5
ptsphy=3600.
pt        = np.zeros((nproma,nlev  ,nblocks), order='F')
pq        = np.zeros((nproma,nlev  ,nblocks), order='F')
pap       = np.zeros((nproma,nlev  ,nblocks), order='F')
paph      = np.zeros((nproma,nlev+1,nblocks), order='F')
plu       = np.zeros((nproma,nlev  ,nblocks), order='F')
plude     = np.zeros((nproma,nlev  ,nblocks), order='F')
pmfu      = np.zeros((nproma,nlev  ,nblocks), order='F')
pmfd      = np.zeros((nproma,nlev  ,nblocks), order='F')
pa        = np.zeros((nproma,nlev  ,nblocks), order='F')
pclv      = np.zeros((nproma,nlev  ,ndim, nblocks), order='F')
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
#dwarf.do_dwarf_call_full(numomp, nproma, nlev, ngptot, ngptotg, nblocks, ptsphy)
print("Filling with 33")
buffer_loc.fill(-33.)
print("Filled with 33")
dwarf.examine_ndarray_flags(buffer_cml)
dwarf.examine_ndarray_flags(buffer_loc)
dwarf.examine_ndarray_flags(pt)
dwarf.examine_ndarray_flags(pt)
rootpath = Path(__file__).resolve().parents[2]
input_path = rootpath/'config-files/input.h5'
input_fields = load_input_fields(path=input_path)
yrecldp, yrmcst, yrethf, yrephli, yrecld = load_input_parameters(path=input_path)

klon = input_fields['KLON']
klev = input_fields['KLEV']

# Get referennce solution fields from file
ref_path = rootpath/'config-files/reference.h5'
ref_fields = load_reference_fields(path=ref_path)

# Populate kernel inputs with raw fields (this splits compound arrays)
satur_args, cloudsc_args = arguments_from_fields(input_fields)
print (cloudsc_args.keys())
#dwarf.do_dwarf_inittest_call(numomp,nproma,nlev,ngptot,nblocks,ngptotg,
#                         ptsphy,
#                         pt,pq,
#                         buffer_cml,buffer_loc,
#                         pap, paph,
#                         plu, plude, pmfu, pmfd,
#                         pa,pclv,psupsat,
#                         pcovptot,
#                         pfplsl, pfplsn, pfhpsl, pfhpsn)
dwarf.do_dwarf_init_call(numomp,nproma,nlev,ngptot,nblocks,ngptotg,
                         ptsphy,
                         pt,pq,
                         buffer_cml,buffer_loc,
                         pap, paph,
                         plu, plude, pmfu, pmfd,
                         pa,pclv,psupsat,
                         pcovptot)

dwarf.do_dwarf_full_call(numomp,nproma,nlev,ngptot,nblocks,ngptotg,
                         ptsphy,
                         pt,pq,
                         buffer_cml,buffer_loc,
                         pap, paph,
                         plu, plude, pmfu, pmfd,
                         pa,pclv,psupsat,
                         pcovptot,
                         pfplsl, pfplsn, pfhpsl, pfhpsn)
dwarf.do_dwarf_validate_call(numomp, nproma, nlev, ngptot, nblocks, ngptotg,
                         ptsphy,
                         pt,pq,
                         buffer_cml,buffer_loc,
                         pap, paph,
                         plu, plude, pmfu, pmfd,
                         pa,pclv,psupsat,
                         pcovptot,
                         pfplsl, pfplsn, pfhpsl, pfhpsn)
