import numpy as np
import freecad.woodWorkbench.tools.WoodMeshGenTools_v11 as WoodMeshGen

def outputParams(nconnector_l,nconnector_l_precrack,nconnector_t,nconnector_t_precrack,\
                 cellwallthickness_early,cellwallthickness_late,height_connector_t,nbeamElem,\
                    long_connector_ratio,segment_length):
    
    # IGA beam parameters
    ninstance = 1
    nsvars_beam = 27 # number of svars for beam
    nsecgp = 4 # number of gauss points for beam sectional integration
    nsvars_secgp = 16 # number of svars at each gauss point
    iprops_beam = np.zeros(6)
    iprops_beam[0]    = 1 # section type index
    iprops_beam[1]    = 1 # number of instances
    iprops_beam[2]    = nconnector_t-nconnector_t_precrack # number of transverse connectors
    iprops_beam[3]    = nconnector_t_precrack # number of precracked transverse connectors
    iprops_beam[4]    = nconnector_l-nconnector_l_precrack # number of longitudinal connectors
    iprops_beam[5]    = nconnector_l_precrack # number of precracked longitudinal connectors
    iprops_beam = [int(x) for x in iprops_beam] 

    props_beam = np.zeros(16)
    props_beam[0]     = 1.5E-7 # Wood substance density [tonne/mm^3]
    props_beam[1]     = 0.8E+4 # Mesoscale elastic modulus [MPa]
    props_beam[2]     = 0.3 # Beam Poisson's ratio
    props_beam[3]     = cellwallthickness_early # Cross-sectional height [mm]
    props_beam[4]     = cellwallthickness_early # Cross-sectional width [mm]
    props_beam[5]     = 100 # Tensile Strength [MPa]
    props_beam[6]     = 200 # Tensile fracture energy [mJ/mm^2]
    props_beam[7]     = 4.1 # Shear Strength Ratio
    props_beam[8]     = 0.2 # Softening Exponent
    props_beam[9]     = 0.2 # Initial Friction
    props_beam[10]    = 0.0 # Asymptotic Friction
    props_beam[11]    = 600 # Transitional Stress [MPa]
    props_beam[12]    = 0.0 # Tensile Unloading
    props_beam[13]    = 0.0 # Shear Unloading
    props_beam[14]    = 0.0 # Shear Softening
    props_beam[15]    = 1.0 # Elastic Analysis Flag

    # Transverse connector parameters
    nsvars_conn_t = 32  # number of svars for transverse connector
    iprops_connector_t_bot = np.zeros(7)
    props_connector_t_bot = np.zeros(26)

    props_connector_t_bot[0]     = 1.5E-7 # Wood substance density [tonne/mm^3]
    props_connector_t_bot[1]     = 0.8E+4 # Mesoscale elastic modulus [MPa]
    props_connector_t_bot[2]     = 0.25 # Shear-Normal coupling coefficient
    props_connector_t_bot[3]     = height_connector_t # Connector height [mm]
    props_connector_t_bot[4]     = 0 # M-Distance [mm]
    props_connector_t_bot[5]     = -height_connector_t/2 # L-Distance [mm]
    props_connector_t_bot[6]     = 30.0 # Tensile Strength [MPa]
    props_connector_t_bot[7]     = 1.0 # Tensile characteristic length [mm] will be updated to # Tensile fracture energy [mJ/mm^2]
    props_connector_t_bot[8]     = 2.6 # Shear Strength Ratio
    props_connector_t_bot[9]     = 0.2 # Softening Exponent
    props_connector_t_bot[10]    = 0.2 # Initial Friction
    props_connector_t_bot[11]    = 0.0 # Asymptotic Friction
    props_connector_t_bot[12]    = 600 # Transitional Stress [MPa]
    props_connector_t_bot[12]    = 0.0 # Tensile Unloading
    props_connector_t_bot[14]    = 0.0 # Shear Unloading
    props_connector_t_bot[15]    = 0.0 # Shear Softening
    props_connector_t_bot[16]    = 0.0 # Elastic Analysis Flag
    props_connector_t_bot[17]    = 0.2 # Compressive Yielding Strength [MPa]
    props_connector_t_bot[18]    = 600 # Initial Hardening Modulus Ratio
    props_connector_t_bot[19]    = 0.0 # Transitional Strain Ratio
    props_connector_t_bot[20]    = 0.0 # Deviatoric Strain Threshold Ratio
    props_connector_t_bot[21]    = 0.0 # Deviatoric Damage Parameter
    props_connector_t_bot[22]    = 0.0 # Final Hardening Modulus Ratio
    props_connector_t_bot[23]    = 0.0 # Densification Ratio
    props_connector_t_bot[24]    = 0.0 # Volumetric Deviatoric Coupling
    props_connector_t_bot[25]    = 0.0 # Compressive Unloading

    iprops_connector_t_bot = [3,ninstance,nbeamElem,nconnector_t-nconnector_t_precrack,nconnector_t_precrack,nconnector_l-nconnector_l_precrack,nconnector_l_precrack]
    iprops_connector_t_bot = [int(x) for x in iprops_connector_t_bot] 
    iprops_connector_t_top = [2,ninstance,nbeamElem,nconnector_t-nconnector_t_precrack,nconnector_t_precrack,nconnector_l-nconnector_l_precrack,nconnector_l_precrack]
    iprops_connector_t_top = [int(x) for x in iprops_connector_t_top] 
    iprops_connector_t_reg = [1,ninstance,nbeamElem,nconnector_t-nconnector_t_precrack,nconnector_t_precrack,nconnector_l-nconnector_l_precrack,nconnector_l_precrack]
    iprops_connector_t_reg = [int(x) for x in iprops_connector_t_reg] 

    props_connector_t_reg = np.copy(props_connector_t_bot)
    props_connector_t_reg[3] = height_connector_t*2
    props_connector_t_reg[5] = 0
    props_connector_t_top = np.copy(props_connector_t_bot)
    props_connector_t_top[5] = height_connector_t/2

    # Longitudinal connector parameters
    nsvars_conn_l = 32  # number of svars for transverse connector
    iprops_connector_l = np.zeros(7)
    props_connector_l = np.zeros(24)

    props_connector_l[0]     = 1.5E-7 # Wood substance density [tonne/mm^3]
    props_connector_l[1]     = 0.8E+4 # Mesoscale elastic modulus [MPa]
    props_connector_l[2]     = 0.18 # Shear-Normal coupling coefficient
    props_connector_l[3]     = cellwallthickness_early*cellwallthickness_late # Connector sectional area [mm^2]
    props_connector_l[4]     = 3.0E+2 # Tensile Strength [MPa]
    props_connector_l[5]     = 0.2E+1# 0.0105 # Tensile characteristic length [mm] will be updated to # Tensile fracture energy [mJ/mm^2]
    props_connector_l[6]     = 2.6 # Shear Strength Ratio
    props_connector_l[7]     = 0.2 # Softening Exponent
    props_connector_l[8]     = 0.2 # Initial Friction
    props_connector_l[9]     = 0.0 # Asymptotic Friction
    props_connector_l[10]    = 600 # Transitional Stress [MPa]
    props_connector_l[11]    = 0.0 # Tensile Unloading
    props_connector_l[12]    = 0.0 # Shear Unloading
    props_connector_l[13]    = 0.0 # Shear Softening
    props_connector_l[14]    = 1.0 # Elastic Analysis Flag
    props_connector_l[15]    = 0.2 # Compressive Yielding Strength [MPa]
    props_connector_l[16]    = 600 # Initial Hardening Modulus Ratio
    props_connector_l[17]    = 0.0 # Transitional Strain Ratio
    props_connector_l[18]    = 0.0 # Deviatoric Strain Threshold Ratio
    props_connector_l[19]    = 0.0 # Deviatoric Damage Parameter
    props_connector_l[20]    = 0.0 # Final Hardening Modulus Ratio
    props_connector_l[21]    = 0.0 # Densification Ratio
    props_connector_l[22]    = 0.0 # Volumetric Deviatoric Coupling
    props_connector_l[23]    = 0.0 # Compressive Unloading

    iprops_connector_l = [1,ninstance,nbeamElem,nconnector_t-nconnector_t_precrack,nconnector_t_precrack,nconnector_l-nconnector_l_precrack,nconnector_l_precrack]
    iprops_connector_l = [int(x) for x in iprops_connector_l] 

    # Strain rate effect parameters
    props_strainrate = np.zeros(4)
    props_strainrate[0] = 0.0 # Strain rate effect flag 
    props_strainrate[1] = 5.0 # Physical time scaling factor
    props_strainrate[2] = 1.0E-5 # Strain rate effect constant c0
    props_strainrate[3] = 5.0E-2 # Strain rate effect constant c1


    return nsvars_beam, nsecgp, nsvars_secgp, iprops_beam, props_beam, \
        nsvars_conn_t, iprops_connector_t_bot, props_connector_t_bot, iprops_connector_t_reg, props_connector_t_reg, \
        iprops_connector_t_top, props_connector_t_top, \
        nsvars_conn_l, iprops_connector_l, props_connector_l, props_strainrate
