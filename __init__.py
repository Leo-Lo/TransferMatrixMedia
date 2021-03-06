from NearFieldOptics.Materials.material_types import *
from NearFieldOptics.Materials.TransferMatrixMedia import MatrixBuilder as mb
from NearFieldOptics.Materials.TransferMatrixMedia import Calculator
# import MatrixBuilder as mb
# import Calculator

#Same private helper method as the one in NearFieldOptics.Material.material_types.
def _prepare_freq_and_q_holder_(freq,q,\
                                angle=None,\
                                entrance=None):
    
    if angle!=None:
        if not entrance: entrance=Air
        angle_rad=angle/180.*pi
        k=safe_sqrt(entrance.optical_constants(freq))*freq
        q=numpy.real(k*numpy.sin(angle_rad))
    else:
        ##Prepare AWA if there are axes in *freq* and *q*##
        freq,q=numerics.broadcast_items(freq,q)
        
    axes=[]; axis_names=[]
    if isinstance(freq,numpy.ndarray):
        freqaxis=freq.squeeze()
        if not freqaxis.ndim: freqaxis.resize((1,))
        axes.append(freqaxis); axis_names.append('Frequency (cm$^{-1}$)')
    if isinstance(q,numpy.ndarray) and angle is None:
        qaxis=q.squeeze()
        if not qaxis.ndim: qaxis.resize((1,))
        axes.append(qaxis); axis_names.append('q-vector (cm$^{-1}$)')
    if axes:
        shape=[len(axis) for axis in axes]
        holder=AWA(numpy.zeros(shape),axes=axes,axis_names=axis_names)
    else: holder=0
    
    return freq,q,ensure_complex(holder)

class LayeredMediaTM(LayeredMedia):
    
    def __init__(self,*layers,layerArrayGUIInput=None,**kwargs):
         
        if layerArrayGUIInput==None:
            self.set_layers(*layers)
         
        else:
            print(layerArrayGUIInput!= None)
            print(layerArrayGUIInput)
            self.set_layers(*layerArrayGUIInput)
        
        #Set default entrance/exit materials
        exkwargs=misc.extract_kwargs(kwargs,entrance=Air,exit=Air)
        self.set_entrance(exkwargs['entrance'])
        self.set_exit(exkwargs['exit'])
         
        self.T_p = mb.TransferMatrix(self,polarization='p')
        self.T_s = mb.TransferMatrix(self,polarization='s')
        
    def reflection_s(self,freq,q=0,angle=None,\
                     entrance=None,exit=None,**kwargs):
        """Get numerical reflection coefficient for s-polarized light.
        
        First the analytical expression for reflection coefficient is assembled. 
        Then the numerical values are evaluated in the helper class Calculator. 
        
        Args:
            freq (array): numpy.ndarray array of frequencies of incident light; in unit of cm^-1
            q (array): numpy.ndarray of in-plane momenta associated with incident light; in unit of cm^-1
        
        Return:
            Numerical reflection coefficient with corresponding dimension of 
            array (based on dimension of freq and q).
        
        """
        freq,q,rsAWA = _prepare_freq_and_q_holder_(freq,q,angle=angle,entrance=entrance)
        C = Calculator.Calculator(self.T_s)
        C.assemble_analytical_reflection_coefficient()
        rs = C.get_numerical_reflection_coefficient(freq,q)
        rsAWA+=rs
        return rsAWA.T
        
    def analytical_reflection_s(self):
        """Get sympy analytical expression of reflection coefficient for s-polarized light."""
        C = Calculator.Calculator(self.T_s)
        C.assemble_analytical_reflection_coefficient()
        rs = C.get_analytical_reflection_coefficient()
        return rs
        
    def reflection_p(self,freq,q=0,angle=None,\
                     entrance=None,exit=None,**kwargs):
        """Get numerical reflection coefficient for p-polarized light.
        
        First the analytical expression for reflection coefficient is assembled. 
        Then the numerical values are evaluated in the helper class Calculator. 
        
        Args:
            freq (array): numpy.ndarray array of frequencies of incident light; in unit of cm^-1
            q (array): numpy.ndarray of in-plane momenta associated with incident light; in unit of cm^-1
        
        Return:
            The numerical reflection coefficient with corresponding dimension of 
            array (based on dimension of freq and q).
        
        """        
        freq,q,rpAWA = _prepare_freq_and_q_holder_(freq,q,angle=angle,entrance=entrance)
        C = Calculator.Calculator(self.T_p)
        C.assemble_analytical_reflection_coefficient()
        rp = C.get_numerical_reflection_coefficient(freq,q)
        rpAWA+=rp
        return rpAWA.T
    
    def analytical_reflection_p(self):
        """Get sympy analytical expression of reflection coefficient for p-polarized light."""
        C = Calculator.Calculator(self.T_p)
        C.assemble_analytical_reflection_coefficient()
        rp = C.get_analytical_reflection_coefficient()
        return rp
    
    def transmission_s(self,freq,q=0,angle=None,\
                     entrance=None,exit=None,**kwargs):
        """Get numerical transmission coefficient for s-polarized light.
        
        First the analytical expression for transmission coefficient is assembled. 
        Then the numerical values are evaluated in the helper class Calculator. 
        
        Args:
            freq (array): numpy.ndarray array of frequencies of incident light; in unit of cm^-1
            q (array): numpy.ndarray of in-plane momenta associated with incident light; in unit of cm^-1
        
        Return:
            Numerical transmission coefficient with corresponding dimension of 
            array (based on dimension of freq and q).
        
        """
        freq,q,tsAWA = _prepare_freq_and_q_holder_(freq,q,angle=angle,entrance=entrance)
        C = Calculator.Calculator(self.T_s)
        C.assemble_analytical_transmission_coefficient()
        ts = C.get_numerical_transmission_coefficient(freq,q)
        tsAWA+=ts
        return tsAWA.T
    
    def analytical_transmission_s(self):
        """Get sympy analytical expression of transmission coefficient for s-polarized light."""
        C = Calculator.Calculator(self.T_s)
        C.assemble_analytical_transmission_coefficient()
        ts = C.get_analytical_transmission_coefficient()
        return ts
        
    def transmission_p(self,freq,q=0,angle=None,\
                     entrance=None,exit=None,**kwargs):
        """Get numerical transmission coefficient for p-polarized light.
        
        First the analytical expression for transmission coefficient is assembled. 
        Then the numerical values are evaluated in the helper class Calculator. 
        
        Args:
            freq (array): numpy.ndarray array of frequencies of incident light; in unit of cm^-1
            q (array): numpy.ndarray of in-plane momenta associated with incident light; in unit of cm^-1
        
        Return:
            The numerical transmission coefficient with corresponding dimension of 
            array (based on dimension of freq and q).
        
        """
        freq,q,tpAWA = _prepare_freq_and_q_holder_(freq,q,angle=angle,entrance=entrance)
        C = Calculator.Calculator(self.T_p)
        C.assemble_analytical_transmission_coefficient()
        tp = C.get_numerical_transmission_coefficient(freq,q)
        tpAWA+=tp
        return tpAWA.T
    
    def analytical_transmission_p(self):
        """Get sympy analytical expression of transmission coefficient for p-polarized light."""
        C = Calculator.Calculator(self.T_p)
        C.assemble_analytical_transmission_coefficient()
        tp = C.get_analytical_transmission_coefficient()
        return tp
    
    def h_field(self,freq,q=0,index=1,angle=None,\
                     entrance=None,exit=None,**kwargs):
        freq,q,hAWA = _prepare_freq_and_q_holder_(freq,q,angle=angle,entrance=entrance)
        C = Calculator.Calculator(self.T_p)
        C.assemble_analytical_H_field(index,'before')
        h = C.get_numerical_H_field(freq,q)
        hAWA+=h
        return hAWA.T
    
    def analytical_h_field(self,index,side):
        C = Calculator.Calculator(self.T_p)
        C.assemble_analytical_H_field(index,side)
        h = C.get_analytical_H_field()
        return h
        
    def Coulomb_kernel(self,freq,q=0,index=1,angle=None,\
                     entrance=None,exit=None,**kwargs):
        freq,q,kAWA = _prepare_freq_and_q_holder_(freq,q,angle=angle,entrance=entrance)
        C = Calculator.Calculator(self.T_p)
        C.assemble_analytical_kernel(index,'before')
        k = C.get_numerical_kernel(freq,q)
        kAWA+=k
        return kAWA.T
        
    def analytical_Coulomb_kernel(self,index,side):
        
        C = Calculator.Calculator(self.T_p)
        C.assemble_analytical_kernel(index,side)
        k = C.get_analytical_kernel()
        return k