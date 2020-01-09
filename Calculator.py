'''
Created on Apr 3, 2019

@author: Leo Lo
'''

from NearFieldOptics.Materials.material_types import *
from NearFieldOptics.Materials.TransferMatrixMedia import MatrixBuilder as mb
import sympy
import copy
import numpy as np
from common.baseclasses import ArrayWithAxes as AWA

class Calculator():
    """Calculator class calculates analytical expression and numerical value of various optical parameters.
        
    Attributes:
        The analytical expression of the following optical parameters are stored:
        - Reflection Coefficient
        - Reflectance
        - Transmission Coefficient
        - Transmittance
        - H Field
        - H field profile
        - E field (x direction) profile
        - E field (z direction) profile
        - Reference Kernel (from Alonso-Gonzalez et al., Nature Nanotechnology 185, 2016)
        - Kernel
        - Number of layers of the LayeredMedium
    """
    
    def __init__(self,transferMatrix):
        """Construct a calculator object. 
        
        Args:
            transferMatrix (TransferMatrix): a transferMatrix object constructed by the MatrixBuilder.py module, based on 
        the input material.
    
        Return:
            void
    
        """
        self.transferMatrix = transferMatrix
        self.analyticalReflectionCoefficient = None
        self.analyticalReflectance = None
        self.analyticalTransmissionCoefficient = None
        self.analyticalTransmittance = None
        self.analyticalHField = None
        self.H_field_profile = None
        self.Ex_field_profile = None
        self.Ez_field_profile = None
        self.analyticalReferenceKernel = None
        self.analyticalKernel = None
        self.numLayers = self.transferMatrix.get_layer_count()-2
    
    def assemble_analytical_reflection_coefficient(self):
        """Create an analytical expression for reflection coefficient of the entire LayeredMedia material.
        
        Args:
            None
        
        Return:
            void
        
        """
        matrix = self.transferMatrix.get_matrix()
        M11 = matrix[0,0]
        M21 = matrix[1,0]
        self.analyticalReflectionCoefficient = M21/M11
    
    def get_analytical_reflection_coefficient(self):
        """Get class variable analyticalReflectionCoefficient.
        
        Args:
             None
        
        Return:
            Analytical expression for reflection coefficient.
        
        """
        return copy.copy(self.analyticalReflectionCoefficient)
    
    def get_numerical_reflection_coefficient(self, freq, q):
        """Get numerical reflection coefficient.
            
        Use lambdify function to substitute numerical values into analytical expression stored in 
        self.analyticalReflectionCoefficient class variable. 
        Broadcast the 1D freq and q arrays into a 2D array to evaluate reflection coefficient at each combination of freq and q.
        
        Args:
            freq (array): numpy.ndarray array of frequencies of incident light; in unit of cm^-1
            q (array): numpy.ndarray of in-plane momenta associated with incident light
        
        Return:
            The numerical reflection coefficient with corresponding dimension of array (based on dimension of freq and q). 
        
        """
        r = self.analyticalReflectionCoefficient
        r_num = self._numerical_evaluation_(r, freq, q)
        return r_num
    
    def assemble_analytical_reflectance(self):
        """Create an analytical expression for reflectance of the entire LayeredMedia material.
        
        Reflectance is the same for both p- and s-polarized lights.
        
        Args:
             None
        
        Return:
            void
        
        """
        self.analyticalReflectance = abs(self.analyticalReflectionCoefficientCoefficient)**2
    
    def get_analytical_reflectance(self):
        """Get class variable analyticalReflectance.
        
        Args:
             None
        
        Return:
            Analytical expression for reflectance.
        
        """
        return copy.copy(self.analyticalReflectance)
    
    def get_numerical_reflectance(self, freq, q):
        """Get numerical reflectance. 
             
        Use lambdify function to substitute numerical values into analytical expression stored in 
        self.analyticalReflectance class variable. 
        
        Args:
            freq (array): numpy.ndarray of frequencies of incident light; in unit of cm^-1
            q (array): numpy.ndarray of in-plane momenta associated with incident light
        
        Return:
            The numerical reflectance with corresponding dimension of array (based on dimension of freq and q).
            
        """
        R = self.analyticalReflectance
        R_num = self._numerical_evaluation_(R, freq, q)
        return R_num
    
    def assemble_analytical_transmission_coefficient(self):
        """Create an analytical expression for transmission coefficient of the entire LayeredMedia material.
        
        Args:
             None
        
        Return:
            void
        
        """
        matrix = self.transferMatrix.get_matrix()
        M11 = matrix[0,0]
        self.analyticalTransmissionCoefficient = 1/M11
    
    def get_analytical_transmission_coefficient(self): 
        """Get class variable analyticalTranmissionCoefficient.
        
        Args:
             None
        
        Return:
            Analytical expression for transmission coefficient.
        
        """
        return copy.copy(self.analyticalTransmissionCoefficient)
    
    def get_numerical_transmission_coefficient(self, freq, q): 
        """Get numerical transmission coefficient. 
             
        Use lambdify function to substitute numerical values into analytical expression stored in 
        self.analyticalTransmissionCoefficient class variable. 
        
        Args:
            freq (array): numpy.ndarray of frequencies of incident light; in unit of cm^-1
            q (array): numpy.ndarray of in-plane momenta associated with incident light
        
        Return:
            The numerical transmission coefficient with corresponding dimension of array (based on dimension of freq and q). 
        
        """
        t = self.analyticalTransmissionCoefficient
        t_num = self._numerical_evaluation_(t, freq, q)
        return t_num
    
    def assemble_analytical_transmittance(self):
        """Create an analytical expression for transmittance of the entire LayeredMedia material.
        
        Based on whether light is p-polarized or s-polarized (info stored in transferMatrix). 
        
        Args:
             None
        
        Return:
            void
        
        """
        epsilon_first,epsilon_last,kz_first,kz_last = sympy.symbols('epsilon_1,epsilon_{0},k_z1,k_z{0}'.format(self.numLayers+2))
        if self.transferMatrix.polarization == 'p':
            self.analyticalTransmittance = epsilon_first*kz_last/(epsilon_last*kz_first)*abs(self.analyticalTransmissionCoefficient)**2
        
        else:       #self.transferMatrix.polarization == 's':
            self.analyticalTransmittance = kz_last/kz_first*abs(self.analyticalTransmissionCoefficient)**2
            
    def get_analytical_transmittance(self):
        """Get class variable analyticalTranmittance.
        
        Args:
             None
        
        Return:
            Analytical expression for transmittance.
        
        """
        return copy.copy(self.analyticalTransmittance)
    
    def get_numerical_transmittance(self, freq, q):
        """Get numerical transmittance. 
             
        Use lambdify function to substitute numerical values into analytical expression stored in 
        self.analyticalTransmittance class variable. 
        
        Args:
            freq (array): numpy.ndarray of frequencies of incident light; in unit of cm^-1
            q (array): numpy.ndarray of in-plane momenta associated with incident light
        
        Return:
            The numerical transmittance with corresponding dimension of array (based on dimension of freq and q). 
        
        """
        T = self.analyticalTransmittance
        T_num = self._numerical_evaluation_(T, freq, q)
        return T_num
    
    def assemble_analytical_H_field(self, n, side):
        """Create analytical expression of H field at either side of the n,n+1 interface; store as a class variable. 
        
        Args:
            n (int): n means that a test charge is placed at the n,n+1 interface. Each layer is indexed; 
                the entrance material has n = 1. Therefore, for a material with N layers, the index goes from 1 to N+2.
            side (str): the side of the n,n+1 interface can be either "before" or "after". The H field
                on the corresponding side is then calculated.
            
        Return:
            void
        
        """
        matrixDictionary = self.transferMatrix.matrixDictionary
        
        #check for parameter inputs
        if n > (self.numLayers+1):
            Logger.raiseException('Index exceed number of layers. n cannot be greater than {0}'.format(self.numLayers+1),exception=ValueError)
        elif n < 1:
            Logger.raiseException('Invalid index. n cannot be less than 1',exception=ValueError)
        elif side!='before' and side!='after':
            Logger.raiseException('The input to side has to either be \'before\' or \'after\'', exception=ValueError)
        
        #begin assembling matrices
        M_1_to_n = sympy.Matrix([[1,0],[0,1]])
        for x in range(2, n+1):
            M_1_to_n *= matrixDictionary["T{0}{1}".format(x-1,x)].get_matrix()
            M_1_to_n *= matrixDictionary["P{0}".format(x)].get_matrix()
            
        M_1_to_n_inv = sympy.Matrix([[M_1_to_n[1,1],-M_1_to_n[0,1]],[-M_1_to_n[1,0],M_1_to_n[0,0]]])
        
        M_n_to_end = mb.TransmissionMatrix(self.transferMatrix.polarization,n,surfaceCurrent='self').get_matrix()
        for x in range(n+1, self.numLayers+2):
            M_n_to_end *= matrixDictionary["P{0}".format(x)].get_matrix()
            M_n_to_end *= matrixDictionary["T{0}{1}".format(x,x+1)].get_matrix()
        
        beta1 = M_1_to_n_inv[0,1]
        delta1 = M_1_to_n_inv[1,1]
        alpha2 = M_n_to_end[0,0]
        gamma2 = M_n_to_end[1,0]
        
        c = sympy.symbols('c')
        
        J = mb.CurrentDensityVector().get_vector()
        inhomogeneousTerm = 4*sympy.pi/c*J/2
        b1 = 1/(beta1*gamma2-alpha2*delta1)*(gamma2*inhomogeneousTerm[0]-alpha2*inhomogeneousTerm[1])
        HfieldBefore = M_1_to_n_inv*sympy.Matrix([[0],[b1]])
        
        if side=='before':
            self.analyticalHField = HfieldBefore
            
        else: #side=='after'
#             transmission = MatrixBuilder.TransmissionMatrix(self.transferMatrix.polarization,n,surfaceCurrent='self').get_matrix()
#             transmission_inv = sympy.Matrix([[transmission[1,1],-transmission[0,1]],[-transmission[1,0],transmission[0,0]]])
#             self.analyticalHField = transmission_inv*(HfieldBefore-inhomogeneousTerm)
            M_nplus1_to_end = sympy.Matrix([[1,0],[0,1]])
            for x in range(n+1,self.numLayers+2):
                M_nplus1_to_end *= matrixDictionary["P{0}".format(x)].get_matrix()
                M_nplus1_to_end *= matrixDictionary["T{0}{1}".format(x,x+1)].get_matrix()
            a_end = 1/(beta1*gamma2-alpha2*delta1)*(delta1*inhomogeneousTerm[0]-beta1*inhomogeneousTerm[1])
            self.analyticalHField = M_nplus1_to_end*sympy.Matrix([[a_end],[0]])
    
    def get_analytical_H_field(self): 
        """Get class variable analyticalHField
        
        Args:
             None.
        
        Return:
            H field right after the (n-1,n) interface.
            
        """
        return copy.copy(self.analyticalHField)
    
    def get_numerical_H_field(self, freq, q): 
        """Get numerical H field. 
            
        Use lambdify function to substitute numerical values into analytical expression stored in 
        self.analyticalHField class variable. 
        
        Args:
            freq (array): numpy.ndarray of frequencies of incident light; in unit of cm^-1
            q (array): numpy.ndarray of in-plane momenta associated with incident light
        
        Return:
            H field with corresponding dimension of array (based on dimension of freq and q). 
        
        """
        H = self.analyticalHField[0] + self.analyticalHField[1]
        H_num = self._numerical_evaluation_(H, freq, q)
        return H_num
    
    def _get_interface_position_list_(self):
        T = self.transferMatrix
        thickness = 0
        list = [0]
        num_layer = T.layerIndex-2
        for i in range(2,num_layer+2):
            d = T.layerDictionary['L'+str(i)].get_thickness()
            thickness += d
            list = np.append(list,thickness)
        return list
    
    def _extract_singleton_array_value_(self,param):
        if isinstance(param,np.ndarray):
            return np.ndarray.item(param)
        return param
    
    def _get_interface_H_field_(self,freq,q,H_0):
        """A helper method to calculate the numerical amplitude of H field right after each interface of the LayeredMedium.
            
        Args:
            freq (array): numpy.ndarray of frequencies of incident light; in unit of cm^-1
            q (array): numpy.ndarray of in-plane momenta associated with incident light
            H_0: a numpy column matrix with two elements, describing the amplitude of H field right before the entrance interface 
            propagating in both directions of the LayeredMedium.
        
        Return:
            amplitude of H field right after each interfaces
        
        """
        
        T = self.transferMatrix
        interface_H_field_array = []
        
        ##entrance interface
        analytical_transmission_matrix = T.matrixDictionary['T12'].get_matrix()
        numerical_transmission_matrix = self._numerical_evaluation_(analytical_transmission_matrix,freq,q)
        tm = np.linalg.inv(numerical_transmission_matrix)
        new_H_field = tm*H_0
        interface_H_field_array.append(new_H_field)
        
        ## rest of the interfaces
        analytical_transfer_matrix = analytical_transmission_matrix
        endIndex = 2+int(len(T.matrixDictionary)/2)
        for i in range(2,endIndex):
            pm = T.matrixDictionary['P'+str(i)].get_matrix()
            tm = T.matrixDictionary['T'+str(i)+str(i+1)].get_matrix()
            analytical_transfer_matrix = analytical_transfer_matrix*pm*tm
            
            numerical_transfer_matrix = self._numerical_evaluation_(analytical_transfer_matrix,freq,q)
            tm = np.linalg.inv(numerical_transfer_matrix)
            new_H_field = tm*H_0
            interface_H_field_array.append(new_H_field)
        
        return interface_H_field_array
    
    def _update_E_field_profile_(self,Ex_profile, Ez_profile, material, H, freq, q):
        
        T = self.transferMatrix
        
        kz = material.get_kz(freq,q)
        epsilon = material.epsilon(freq,q)
        kz = self._extract_singleton_array_value_(kz)
        epsilon = self._extract_singleton_array_value_(epsilon)
        
        omega = 2*np.pi*freq
        Ex = (H[1]-H[0])*29979245368*kz/(omega*epsilon)
        Ex_profile = np.append(Ex_profile,Ex)
        Ez = (H[0]+H[1])*29979245368*q/(omega*epsilon)
        Ez_profile = np.append(Ez_profile,Ez)
        
        return Ex_profile,Ez_profile
    
    def _numerical_evaluation_(self,analytical_quantity,freq,q):
        """Substitute numerical values into any analytical expression.
        
        Use lambdify function to substitute numerical values into analytical quantity
        specified by user.
        Automatically broadcast the 1D freq and q arrays into a 2D array to evaluate reflection coefficient at each combination of freq and q.
        
        Args:
            freq (array): numpy.ndarray array of frequencies of incident light; in unit of cm^-1
            q (array): numpy.ndarray of in-plane momenta associated with incident light
        
        Return:
            The numerical value of analytical_quantity with corresponding dimension of array (based on dimension of freq and q). 
        
        """
        T = self.transferMatrix
        entranceMaterial = T.entrance
        exitMaterial = T.exit
        layerDictionary = T.layerDictionary
    
        subs = {}
        subs['c'] = 3e10
        subs['omega'] = 2*np.pi*freq
    
        #for first boundary
        subs['k_z1'] = entranceMaterial.get_kz(freq,q)
        subs['epsilon_1'] = entranceMaterial.epsilon(freq,q)
        subs['mu_1'] = entranceMaterial.mu(freq,q)
    
        for x in range(2, self.numLayers+2):
    
            layer = layerDictionary['L'+str(x)]
            material = layer.get_material()
            surface = layerDictionary['S'+str(x-1)+str(x)]
            subs['k_z{}'.format(x)] = material.get_kz(freq,q)
            subs['z{}'.format(x)] = layer.get_thickness()
            subs['sigma{0}{1}'.format(x-1,x)] = surface.conductivity(freq)
            subs['epsilon_{}'.format(x)] = material.epsilon(freq,q)
            subs['mu_{}'.format(x)] = material.mu(freq,q)
    
        #for last boundary
        subs['k_z{}'.format(self.numLayers+2)] = exitMaterial.get_kz(freq,q)
        subs['epsilon_{}'.format(self.numLayers+2)] = exitMaterial.epsilon(freq,q)
        subs['mu_{}'.format(self.numLayers+2)] = exitMaterial.mu(freq,q)
        surface = layerDictionary['S'+str(self.numLayers+1)+str(self.numLayers+2)]
        subs['sigma{0}{1}'.format(self.numLayers+1,self.numLayers+2)] = surface.conductivity(freq)
    
        numerics = sympy.lambdify(subs.keys(), analytical_quantity, modules='numpy')
        numerical_quantity = numerics(*subs.values())
        return numerical_quantity
    
    def compute_field_profile(self,freq,q,a=1.,distance_into_entrance=0,distance_into_exit=0,num_sample=1000,subtract_incident_field=True,normalized=True):
        """Calculate the numerical values of E field (z direction) profile, E field (x direction) profile, and H field (y direction) profile. 
        
        Field profile means the value of the field as a function of z position, with z axis normal to interface.
        
        Note: in order to propagate from entrance to exit, due to the formulation used in this code, we need to use
        the inverse of all the matrices.
        
        Args:
            freq (float): Frequency of incident light; in unit of cm^-1
            q (float): In-plane momentum of incident light; in unit of cm^-1
            a (float): magnitude of incident H field
            distance_into_entrance (float): distance of profile before entrance; in unit of cm
            distance_into_exit (float): distance of profile after exit; in unit of cm
            num_sample: number of position to sample fields
        
        """
        #Compute parameters needed for calculating fields
        self.assemble_analytical_reflection_coefficient()
        b=self.get_numerical_reflection_coefficient(freq,q)*a
        T = self.transferMatrix
        num_layer = T.layerIndex-2
        index = 2
        H_0 = np.matrix([[a],[b]])
        H_profile = []
        Ex_profile = []
        Ez_profile = []
        omega = 2*np.pi*freq
        
        #Compute parameters related to position of interfaces
        interface_position_list = self._get_interface_position_list_()
        thickness = interface_position_list[-1]+distance_into_entrance+distance_into_exit
        interface_index = 1    # the index of the next interface
        if len(interface_position_list)>1:
            next_interface_position = interface_position_list[interface_index]
        step_size = thickness/(num_sample-1)
        positionArray = np.linspace(-distance_into_entrance,interface_position_list[-1]+distance_into_exit,num=int(num_sample))
        
        startingIndex = int(distance_into_entrance/step_size)+1
        endingIndex = int(interface_position_list[-1]/step_size)+startingIndex
        
        #Obtain fields inside entrance
        if subtract_incident_field == True:
            H_0 = np.matrix([[0],[b]])
        for z in positionArray[0:startingIndex]:
            material = T.entrance
            kz = material.get_kz(freq,q)
            propagation_matrix = np.matrix([[np.exp(-1j*kz*abs(z)) , 0],
                                            [0 , np.exp(1j*kz*abs(z))]])
            H = propagation_matrix*H_0
            H_profile = np.append(H_profile,H.sum())
            Ex_profile, Ez_profile = self._update_E_field_profile_(Ex_profile, Ez_profile, material, H, freq, q)
        
        #Obtain interface fields
        H_0 = np.matrix([[a],[b]])
        interface_H_field_array = self._get_interface_H_field_(freq,q,H_0)
        H_at_interface = interface_H_field_array[0]
        
        #Obtain fields inside LayeredMedium
        distance_from_interface = step_size
        for z in positionArray[startingIndex:endingIndex]:
        
            floating_error = step_size/1000
            if z-floating_error > next_interface_position:
                #go to next interface
                H_at_interface = interface_H_field_array[interface_index]
                interface_index += 1
                distance_from_interface = z-next_interface_position
                next_interface_position = interface_position_list[interface_index]
        
            material = T.layerDictionary['L'+str(interface_index+1)].get_material()
            kz = material.get_kz(freq,q)
            propagation_matrix = np.matrix([[np.exp(1j*kz*distance_from_interface) , 0],
                                            [0 , np.exp(-1j*kz*distance_from_interface)]])
        
            H = propagation_matrix*H_at_interface
            H_profile = np.append(H_profile,H.sum())
            Ex_profile, Ez_profile = self._update_E_field_profile_(Ex_profile, Ez_profile, material, H, freq, q)
        
            distance_from_interface += step_size
        
        # Obtain fields inside exit
        for z in positionArray[endingIndex:]:
                material = T.exit
                kz = material.get_kz(freq,q)
                propagation_matrix = np.matrix([[np.exp(1j*kz*abs(z)) , 0],
                                                [0 , np.exp(-1j*kz*abs(z))]])
                H_exit = interface_H_field_array[-1]
                H = propagation_matrix*H_exit
                H_profile = np.append(H_profile,H.sum())
                Ex_profile, Ez_profile = self._update_E_field_profile_(Ex_profile, Ez_profile, material, H, freq, q)
        
        #Subtract incident field
        if subtract_incident_field == True:
            H_incident_profile = np.zeros(startingIndex)        #The incident field is already removed when propagating into entrance
            Ex_incident_profile = np.zeros(startingIndex)        #The incident field is already removed when propagating into entrance
            Ez_incident_profile = np.zeros(startingIndex)        #The incident field is already removed when propagating into entrance
            for z in positionArray[startingIndex:]:
                propagation_matrix = np.matrix([[np.exp(1j*kz*abs(z)) , 0],
                                                [0 , np.exp(-1j*kz*abs(z))]])
                H_0_subtracted = np.matrix([[a],[0]])
                H = propagation_matrix * H_0_subtracted
                H_incident_profile = np.append(H_incident_profile,H.sum())
                Ex_incident_profile, Ez_incident_profile = self._update_E_field_profile_(Ex_incident_profile, Ez_incident_profile, material, H, freq, q)
            H_profile = H_profile - H_incident_profile
            Ex_profile = Ex_profile - Ex_incident_profile
            Ez_profile = Ez_profile - Ez_incident_profile
            
            H0 = b
        
        elif subtract_incident_field == False:
            H0=a+b
        else:
            Logger.raiseException('Invalid input for argument \'subtract_incident_field\'. Can only be boolean value.', exception=ValueError)
        
        if normalized==True:
            material = T.entrance
            kz = material.get_kz(freq,q)
            epsilon = material.epsilon(freq,q)
            Ez0 = H0*29979245368*q/(omega*epsilon)
            Ez_profile = Ez_profile/Ez0*b
            Ex0 = H0*29979245368*kz/(omega*epsilon)
            Ex_profile = Ex_profile/Ex0*b
        elif (normalized!=False)&(normalized!=True):
            Logger.raiseException('Invalid input for argument \'normalized\'. Can only be boolean value.', exception=ValueError)
        
        self.Ez_field_profile = AWA(Ez_profile,axes=[positionArray*1e7],axis_names=['distance from entrance (nm)'])
        self.Ex_field_profile = AWA(Ex_profile,axes=[positionArray*1e7],axis_names=['distance from entrance (nm)'])
        self.H_field_profile = AWA(H_profile,axes=[positionArray*1e7],axis_names=['distance from entrance (nm)'])

    def get_H_field_profile(self):
        """Get H field profile. Use after the compute_field_profile method to obtain nonempty result.
        
        Field profile means the value of the field as a function of z position, with z axis normal to interface.
        
        Args:
             None.
        
        Return:
            class variable H_field_profile 
            
        """
        return copy.copy(self.H_field_profile)
    
    def get_Ez_field_profile(self):
        """Get E field (z direction) profile. Use after the compute_field_profile method to obtain nonempty result.
        
        Field profile means the value of the field as a function of z position, with z axis normal to interface.
        
        Args:
             None.
        
        Return:
            class variable Ez_field_profile 
            
        """
        return copy.copy(self.Ez_field_profile)
    
    def get_Ex_field_profile(self):
        """Get E field (x direction) profile. Use after the compute_field_profile method to obtain nonempty result.
        
        Field profile means the value of the field as a function of z position, with z axis normal to interface.
        
        Args:
             None.
        
        Return:
            class variable Ex_field_profile 
            
        """
        return copy.copy(self.Ex_field_profile)
    
    def get_2d_field_profile(self,q,field_str='Ez',num_sample=10000,x_window_size = 4):
        if field_str=='Ez':
            field = self.Ez_field_profile
        elif field_str=='Ex':
            field = self.Ex_field_profile
        elif field_str=='H':
            field = self.H_field_profile
        else:
            Logger.raiseException('Invalid input for field_str argument. Only takes Ez,Ex, or H as input.',exception=ValueError)    
        
        wavelength = 2*np.pi/q
        step_size = x_window_size*wavelength/num_sample
        x_array = np.linspace(0,x_window_size*wavelength,num_sample)
        propagate_array = np.cos(q*x_array)
        field_2d = np.outer(field,propagate_array)  
        return field_2d
    
    
    def assemble_analytical_reference_kernel(self):
        """Create an analytical expression for Coulomb kernel from Alonso-Gonzalez et al., Nature Nanotechnology 185, 2016.
        
        The material is a graphene sheet encapsulated by two uniaxial layers in an isotropic medium (permittivity of 
        medium above and below the material can be different). 
        
        Args:
             None
        
        Return:
            void
        
        """
        epsilon_x,epsilon_z,epsilon_a,epsilon_b,e,q,d1,d2 = sympy.symbols('epsilon_x,epsilon_z,epsilon_a,epsilon_b,e,q,d1,d2')
        v_q = 4*sympy.pi*e**2/(q*(epsilon_a+epsilon_b))
        epsilon_tilta = (epsilon_a*epsilon_b+epsilon_x*epsilon_z)/(epsilon_a+epsilon_b)
        V = v_q*sympy.Rational(1,2)*(sympy.sqrt(epsilon_x*epsilon_z)
                                     + (epsilon_a+epsilon_b)*sympy.tanh(q*sympy.sqrt(epsilon_x/epsilon_z)*(d1+d2))
                                     + (epsilon_b-epsilon_a)*sympy.sinh(q*sympy.sqrt(epsilon_x/epsilon_z)*(d1-d2))/sympy.cosh(q*sympy.sqrt(epsilon_x/epsilon_z)*(d1+d2))
                                     + (sympy.sqrt(epsilon_x*epsilon_z)-epsilon_a*epsilon_b/sympy.sqrt(epsilon_x*epsilon_z))*sympy.cosh(q*sympy.sqrt(epsilon_x/epsilon_z)*(d2-d1))/sympy.cosh(q*sympy.sqrt(epsilon_x/epsilon_z)*(d1+d2))
                                     + epsilon_a*epsilon_b/sympy.sqrt(epsilon_x*epsilon_z)
                                     )/(
                                         sympy.sqrt(epsilon_x*epsilon_z)+epsilon_tilta*sympy.tanh(q*sympy.sqrt(epsilon_x/epsilon_z)*(d1+d2))
                                         )
        
        self.analyticalReferenceKernel = V
    
    def get_analytical_reference_kernel(self): 
        """Get analytical Coulomb kernel from Alonso-Gonzalez et al., Nature Nanotechnology 185, 2016.
        
        The material is a graphene sheet encapsulated by two uniaxial layers in an isotropic medium (permittivity of 
        medium above and below the material can be different). 
            
        Args:
            None.
        
        Return:
            Analytical expression of Coulomb kernel for an graphene encapsulated by two uniaxial materials
            in an isotropic medium.
        
        """
        return copy.copy(self.analyticalReferenceKernel)
    
    def get_numerical_reference_kernel(self,freq,q,material,d1,d2,epsilon_a=1,epsilon_b=1):
        """Get numerical Coulomb kernel from Alonso-Gonzalez et al., Nature Nanotechnology 185, 2016.
            The material is a graphene sheet encapsulated by two uniaxial layers in an isotropic medium.
        
        Args:
            q (float array): an array of in-plane momenta of incident light.
            epsilon_x (float): the complex in-plane relative permittivity of uniaxial material
            epsilon_z (float): the complex out-of-plane relative permittivity of uniaxial material
            epsilon_a (float): the complex relative permittivity of isotropic medium above the sample 
            epsilon_b (float): the complex relative permittivity of isotropic medium below the sample
        
        Return:
            An array of numerical value of Coulomb kernel (as a function of q) for an graphene encapsulated by two 
            uniaxial materials in an isotropic medium.
        
        """
        V = self.analyticalReferenceKernel
        
        subs = {}
        subs['epsilon_a'] = epsilon_a
        subs['epsilon_b'] = epsilon_b
        subs['e'] = 1        #"normalized"
        subs['d1'] = d1
        subs['d2'] = d2
        subs['q'] = q
        
        if (type(material)==BaseIsotropicMaterial or type(material)==IsotropicMaterial):
            subs[sympy.symbols('epsilon_x')] = material.epsilon(freq,q)
            subs[sympy.symbols('epsilon_z')] = material.epsilon(freq,q)
        elif (type(material)==BaseAnisotropicMaterial or type(material)==AnisotropicMaterial):
            subs[sympy.symbols('epsilon_x')] = material.ordinary_epsilon(freq,q)
            subs[sympy.symbols('epsilon_z')] = material.extraordinary_epsilon(freq,q)
        else:
            Logger.raiseException('Invalid material. Accept only material of type BaseIsotropicMaterial,\
            IsotropicMaterial,BaseAnisotropicMaterial, or AnisotropicMaterial.',exception=ValueError)
        
        numerics = sympy.lambdify(subs.keys(), V, modules='numpy')
        potentialArray = numerics(*subs.values())
          
        return potentialArray
    
    def assemble_analytical_reference_kernel_2(self):
        epsilon_x,epsilon_z,epsilon_a,epsilon_b,e,q,d = sympy.symbols('epsilon_x,epsilon_z,epsilon_a,epsilon_b,e,q,d')
        v_q = 4*sympy.pi*e**2/(q*(epsilon_a+epsilon_b))
        V = v_q*sympy.Rational(1,2)*(
                sympy.sqrt(epsilon_x*epsilon_z)+(epsilon_a+epsilon_b)*sympy.tanh(q*d*sympy.sqrt(epsilon_x/epsilon_z))
            )/(
                sympy.sqrt(epsilon_x*epsilon_z)+(epsilon_x*epsilon_z+epsilon_b*epsilon_a)*sympy.tanh(q*d*sympy.sqrt(epsilon_x/epsilon_z))/(epsilon_a+epsilon_b)
                )
        
        self.analyticalReferenceKernel = V
            
    def get_numerical_reference_kernel_2(self,freq,q,material,d,epsilon_a=1,epsilon_b=1):
        
        V = self.analyticalReferenceKernel
        
        subs = {}
        subs['epsilon_a'] = epsilon_a
        subs['epsilon_b'] = epsilon_b
        subs['e'] = 1        #"normalized"
        subs['d'] = d
        subs['q'] = q
        
        if (type(material)==BaseIsotropicMaterial or type(material)==IsotropicMaterial):
            subs['epsilon_x'] = material.epsilon(freq,q)
            subs['epsilon_z'] = material.epsilon(freq,q)
        elif (type(material)==BaseAnisotropicMaterial or type(material)==AnisotropicMaterial):
            subs['epsilon_x'] = material.ordinary_epsilon(freq,q)
            subs['epsilon_z'] = material.extraordinary_epsilon(freq,q)
        else:
            Logger.raiseException('Invalid material. Accept only material of type BaseIsotropicMaterial,\
            IsotropicMaterial,BaseAnisotropicMaterial, or AnisotropicMaterial.',exception=ValueError)
        
        numerics = sympy.lambdify(subs.keys(), V, modules='numpy')
        potentialArray = numerics(*subs.values())
          
        return potentialArray
    
    def direct_numerical_reference_kernel_2(self,freq,q,material,d,epsilon_a=1,epsilon_b=1):
        
        if (type(material)==BaseIsotropicMaterial or type(material)==IsotropicMaterial):
            epsilon_x = material.epsilon(freq,q)
            epsilon_z = material.epsilon(freq,q)
        elif (type(material)==BaseAnisotropicMaterial or type(material)==AnisotropicMaterial):
            epsilon_x = material.ordinary_epsilon(freq,q)
            epsilon_z = material.extraordinary_epsilon(freq,q)
        else:
            Logger.raiseException('Invalid material. Accept only material of type BaseIsotropicMaterial,\
            IsotropicMaterial,BaseAnisotropicMaterial, or AnisotropicMaterial.',exception=ValueError)
        
        e = 1
        v_q = -4*np.pi*e/(q*(epsilon_a+epsilon_b))
        V = v_q*(
                safe_sqrt(epsilon_x*epsilon_z)+epsilon_b*np.tanh(q*d*safe_sqrt(epsilon_x/epsilon_z))
            )/(
                safe_sqrt(epsilon_x*epsilon_z)+(epsilon_x*epsilon_z+epsilon_b*epsilon_a)*np.tanh(q*d*safe_sqrt(epsilon_x/epsilon_z))/(epsilon_a+epsilon_b)
                )
        return V
    
    def assemble_analytical_kernel(self,n,side):
        """Create analytical expression of Coulomb kernel from transfer matrix method.
        
        Position of the kernel is at either side of the n,n+1 interface.
        Analytical kernel is stored as a class variable.
        
        Args:
            n (int): n means that a test charge is placed at the n,n+1 interface. Each layer is indexed; 
                the entrance material has n = 1. Therefore, for a material with N layers, the index goes from 1 to N+2.
            side (str): the side of the n,n+1 interface can be either "before" or "after". The H field
                on the corresponding side is then calculated.
        
        Return:
            void
        
        """
        if side == 'before':
            epsilon_n = sympy.symbols('epsilon_{}'.format(n))
            k_n = sympy.symbols('k_z{}'.format(n))
        
        elif side == 'after':
            epsilon_n = sympy.symbols('epsilon_{}'.format(n+1))
            k_n = sympy.symbols('k_z{}'.format(n+1))
            
        else:
            Logger.raiseException('The input to side has to either be \'before\' or \'after\'', exception=ValueError)
        
        self.assemble_analytical_H_field(n,side)
        omega,c,q = sympy.symbols('omega,c,q')
        a = self.get_analytical_H_field()[0]
        b = self.get_analytical_H_field()[1]
        
        self.analyticalKernel = -sympy.I*c*k_n/(omega*epsilon_n*q)*(b-a)
    
    def get_analytical_kernel(self): 
        """Get analytical Coulomb kernel from transfer matrix method.
        
        Args:
             None
        
        Return:
            analytical expression of the Coulomb kernel (self.analyticalKernel) at the n-1,n interface.
            
        """
        return copy.copy(self.analyticalKernel)
    
    def get_numerical_kernel(self, freq, q):
        """Get numerical Coulomb kernel from transfer matrix method.
            
        Use lambdify function to substitute numerical values into analytical expression stored in 
        self.analyticalHField class variable.
        
        Args:
            freq (array): numpy.ndarray of frequencies of incident light; in unit of cm^-1
            q (array): numpy.ndarray of in-plane momenta associated with incident light
        
        Return:
            The numerical Coulomb kernel with corresponding dimension of array (based on dimension of freq and q).
        
        """

        V = self.analyticalKernel
        V_num = self._numerical_evaluation_(V, freq, q)
        return V_num
        
