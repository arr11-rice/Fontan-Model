import scipy.optimize
import numpy as np

def inverse_compliance(variables, *param):
    """
    Given a desired systemic arterial pressure (P_sa), solve for the compliance values (C_d, C_s, C_sa, C_pv, C_pa).
    """
    (P_sa_target, UVR, LVR, PVR, HR) = param  # Fixed parameters
    (C_d, C_s, C_sa, C_pv, C_pa) = variables  # Unknown compliance values
    
    # Solve for flows using current compliance estimates
    def fun_flows(variables, *param):
        (UVR, LVR, PVR, HR, C_d, C_s, C_sa, C_pv, C_pa) = param
        (Q_v, Q_u, Q_l, Q_p, P_sa, P_pa, P_pv) = variables
        
        eqn_a01 = Q_v - HR * (C_d * P_pv - C_s * P_sa)
        eqn_a02 = Q_u + Q_l - Q_v
        eqn_a03 = Q_p - Q_v
        eqn_a04 = PVR * Q_p - (P_pa - P_pv)
        eqn_a05 = UVR * Q_u - (P_sa - P_pa)
        eqn_a06 = LVR * Q_l - (P_sa - P_pa)
        eqn_a07 = 1 - (C_sa * P_sa + C_pv * P_pv + C_pa * P_pa)
        
        return [eqn_a01, eqn_a02, eqn_a03, eqn_a04, eqn_a05, eqn_a06, eqn_a07]
    
    # Initial guesses for flows
    z0_flows = (3.1, 1.5, 1.5, 3.2, 75, 26, 2)
    
    # Solve for flows given compliance values
    result_flows = scipy.optimize.fsolve(fun_flows, z0_flows, args=(UVR, LVR, PVR, HR, C_d, C_s, C_sa, C_pv, C_pa))
    (Q_v, Q_u, Q_l, Q_p, P_sa, P_pa, P_pv) = result_flows
    
    # Equations to enforce desired P_sa while maintaining consistency across compliance values
    eqn_c1 = P_sa - P_sa_target
    eqn_c2 = C_d * P_pv - C_s * P_sa - (Q_v / HR)
    eqn_c3 = C_sa * P_sa + C_pv * P_pv + C_pa * P_pa - 1
    eqn_c4 = Q_p - Q_v
    eqn_c5 = PVR * Q_p - (P_pa - P_pv)
    
    return [eqn_c1, eqn_c2, eqn_c3, eqn_c4, eqn_c5]

# Define target P_sa and fixed parameters
P_sa_target = 75  # Desired systemic arterial pressure
UVR, LVR, PVR, HR = 60, 40, 10, 150  # Fixed parameters

# Initial guesses for compliance values
z0_compliances = (0.02, 0.0001, 1/135, 30 * (1/135), 2 * (1/135))

# Solve for compliance values that achieve the target P_sa
result_compliances = scipy.optimize.fsolve(inverse_compliance, z0_compliances, args=(P_sa_target, UVR, LVR, PVR, HR))
(C_d, C_s, C_sa, C_pv, C_pa) = result_compliances

print("Optimal compliance values:")
print(f"C_d: {C_d}, C_s: {C_s}, C_sa: {C_sa}, C_pv: {C_pv}, C_pa: {C_pa}")


#testing correctness
def fun_flows_after(variables, *param):
    (UVR, LVR, PVR, HR, C_d, C_s, C_sa, C_pv, C_pa) = param
    (Q_v, Q_u, Q_l, Q_p, P_sa, P_pa, P_pv) = variables
    
    eqn_a01 = Q_v - HR * (C_d * P_pv - C_s * P_sa)
    eqn_a02 = Q_u + Q_l - Q_v
    eqn_a03 = Q_p - Q_v
    eqn_a04 = PVR * Q_p - (P_pa - P_pv)
    eqn_a05 = UVR * Q_u - (P_sa - P_pa)
    eqn_a06 = LVR * Q_l - (P_sa - P_pa)
    eqn_a07 = 1 - (C_sa * P_sa + C_pv * P_pv + C_pa * P_pa)
    
    return [eqn_a01, eqn_a02, eqn_a03, eqn_a04, eqn_a05, eqn_a06, eqn_a07]

param_flows = (UVR, LVR, PVR, HR, C_d, C_s, C_sa, C_pv, C_pa)
z0_flows = (3.1, 1.5, 1.5, 3.2, 75, 26, 2)

result_flows = scipy.optimize.fsolve(fun_flows_after, z0_flows, args=param_flows, full_output=True, xtol=1e-4)
(Q_v, Q_u, Q_l, Q_p, P_sa, P_pa, P_pv) = result_flows[0]
print(P_sa) #should be close to the ideal we defined earlier