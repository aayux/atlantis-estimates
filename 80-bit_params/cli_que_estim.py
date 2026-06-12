###
# ###
# Accompanying script for the paper
# "Practical, Round-Optimal Lattice-Based Blind Signatures "
# The script provides security estimates for the blind signature construction
# as well as estimate for the signature size
# ###
# Security and size estimates derived from [LNP22] ***version 1***
#(dated March 22, 2022), available at https://eprint.iacr.org/2022/
#
####DS: the above needs to be updated
###

from MSIS_security import SIS_optimize_attack, SIS_l2_cost, SIS_linf_cost
from MLWE_security import MLWE_optimize_attack, LWE_primal_cost, LWE_dual_cost
from model_BKZ import svp_classical
from math import sqrt, floor, ceil, log, exp

pi = 3.1415926535897932384626433832795028841971693993751058209749


Kilo = 1024
G_entropy_const = 4.13
#4.13 as in 2022/141

G_tail_const = 2.04614178164472 # 2.57 is in LNP
# this is log(4.13,2) as in 2022/141

###############################################################################
#
# variables that define LWE/SIS hardness differ in [LNP22] from Fig.10 to Fig.6
# below is the transition dictionary.
# We instantiate Fig.10, monitor how the parameters get updated to Fig.6
# and instantiate Thm. 4.2. with the updated parameters. The result is given
# in the HARDNESS section below.
#
#               |  Fig.10    |    Fig.8                 | Fig.6
# ______________________________________________________________________________
#   A1.nrows    |   n        |    n                     |    n
#   s1.dim      |   m1       |    m1+ve                 |  m1+ve
#   s2.dim      |   m2       |    m2                    |    m2
#   m.dim       |   ell      |    ell+(256/d+1)*{0,1,2} | Fig.8+lamb/2
# ______________________________________________________________________________
#                               HARDNESS (Thm. 4.2)
#                       here m1 and ell are those of Fig.10
#   LWE:        |
#      dim      |   m2-n-ell-{0,1,2}*(256/d+1)-lamb/2-1
#      Nsamples |   n+ell+{0,1,2}*(256/d+1)+lamb/2+1
# ______________________________________________________________________________
#   SIS:
#    hight of A |   n
#    width of A |   m1+m2+ve
#
#
#   multiplication by a set in m.dim for Fig.8 and hardness must be
#   understood as follows:  we multiply
#       by 0 if ve=vd=0,
#       by 1 if only one of {ve,vd} is at least 1,
#       by 2 if both {ve,vd} are at least 1
#
#
############################################################################

def SIS_security(paramset):
    # as per [LNP22] Thm 4.2 (see above dictionary)
    d  = paramset['d']
    q  = paramset['q']
    n  = paramset['n']
    m1 = paramset['m1']
    m2 = paramset['m2']
    norm_s1  = paramset['norm_s1']
    ve = paramset['ve']
    nu = paramset['nu']
    eta    = paramset['eta']
    gamma1 = paramset['gamma1']
    gamma2 = paramset['gamma2']
    gamma  = paramset['gamma'] #Dilithium-G - compression
    D  = paramset['D'] #Dilithium-G - compression
    sigma1 = gamma1*eta*sqrt(norm_s1+ve*d)
    # sqrt(norm_s1/2+ve*d)is an upper bound on the ell2-norm of ABDLOP's s1.
    sigma2 = gamma2*eta*nu*sqrt(m2*d)
    # sqrt(m2*nu*d) is an upper bound on the ell2-norm of ABDLOP's s2.

    # Bound on the SIS solution. From LNP22
    B1 = 2*sigma1*sqrt(2*m1*d)
    B2 = 2*sigma2*sqrt(2*m2*d)+2**D*eta*sqrt(n*d)+gamma*sqrt(n*d)
    BSIS = 4*eta*sqrt(B1**2+B2**2)

    assert(BSIS<q)
    # maximal width
    max_w = (m1+m2)*d
    h = n*d
    ##
    # estimation for the SIS problem A*x = 0 mod q, where A is from Zq^(h X max_w) and ||x||_2 < BSIS
    # returns optimal m, beta, bit complexity
    m_pc, b_pc, c_pc = SIS_optimize_attack(q, max_w, h, BSIS, cost_attack=SIS_l2_cost, cost_svp=svp_classical, verbose=False)
    return (m_pc, b_pc, c_pc)


def LWE_security(paramset, attack=LWE_primal_cost):
    # as per [LNP22] Thm 4.2 (see above dictionary)
    d  = paramset['d']
    q  = paramset['q']
    n  = paramset['n']
    m1 = paramset['m1']
    m2 = paramset['m2']
    nu = paramset['nu']
    ve = paramset['ve']
    vd = paramset['vd']
    ell = paramset['ell']
    lamb = paramset['lambda']
    assert(lamb%2 == 0)
    if vd==0 and ve==0:
        scalar = 0
    elif vd>=1 and ve>=1:
        scalar = 2
    else: scalar = 1
    ell_updated = ell+(round(256/d)+1)*scalar+round(lamb/2)
    mLWE = (n+ell_updated+1)*d
    nLWE = (m2-n-ell_updated-1)*d
    ##
    # estimation for the LWE problem
    # returns optimal m, beta, bit complexity
    m_pc_, b_pc_, c_pc_ = MLWE_optimize_attack(q, nLWE, mLWE, nu, cost_attack=attack, cost_svp=svp_classical, verbose=False)
    return(m_pc_, b_pc_, c_pc_)


def proof_size(paramset):
    d  = paramset['d']
    q  = paramset['q']
    n  = paramset['n']
    m1 = paramset['m1']
    m2 = paramset['m2']
    assert(m1>=512/d)
    assert(m2>=512/d)      # as per Thm 4.3
    alphae  = paramset['alphae']
    norm_s1 = paramset['norm_s1']
    ell = paramset['ell']
    lamb   = paramset['lambda']
    assert(lamb%2 == 0)    # as per Section 4.4
    ve = paramset['ve']
    vd = paramset['vd']
    if vd==0 and ve==0:
        scalar = 0
    elif vd>=1 and ve>=1:
        scalar = 2
    else: scalar = 1
    #ell_updated = ell+(round(256/d)+1)*scalar+lamb/2
    #print('ell_updated:',ell_updated)
    ##
    D     = paramset['D']
    nu = paramset['nu']
    eta    = paramset['eta']
    gamma1 = paramset['gamma1']
    gamma2 = paramset['gamma2']
    gammae = paramset['gammae']
    sigma1 = gamma1*eta*sqrt(norm_s1+ve*d) # norm_s1 is the squared norm of ABDLOP's s1 as per Thm. 5.3 (Fig. 10), hence add ve*d (norm of x)
    sigma2 = gamma2*eta*nu*sqrt(m2*d)

    sigmae = gammae*sqrt(337)*(sqrt(alphae))
    # See [LNP22] fig 10; this is an upper bound on the norm of s^(e). Note that the bound in "public information" is flawed (misplaced square roots)
    print('proof-size log-sigmas:', log(sigma1, 2), log(sigma2, 2), log(sigmae, 2))
    lgq             = (log(q,2))
    challnge_size   = ceil(log(2*kappa+1,2))*d
    hint_size = 2.25*n*d
    # p.49 of LNP22 above the paragraph "Dilithium compression" + we already have ve in m1 + we do not send t^{(d)}'s
    size_plain      = (n+ell+(256/d+1)+1+lamb)*d*lgq + m1*d*(G_tail_const+(log(sigma1, 2))) + m2*d*(G_tail_const+(log(sigma2,2)))+ 256*(G_tail_const+(log(sigmae,2))) + challnge_size
    # p.50 (top) of LNP22
    size_cut        = n*d*(lgq - D)+(ell+(256/d+1)+1+lamb)*d*lgq + m1*d*(G_tail_const+(log(sigma1, 2))) + (m2-n)*d*(G_tail_const+(log(sigma2,2)))+ 256*(G_tail_const+(log(sigmae,2))) + challnge_size + hint_size

    """ OLD CODE """
    ##
    # instead of sending full t_A = t_A0+2^D*t_A1, we only send the high order bits t_A1
    # As in Dilithium (https://eprint.iacr.org/2017/633.pdf) we send the carry bits that appear by not adding c*t_A0
    # The prover runs MakeHint(-c*t_A0, w+c*t_A0, |c*t_A0|_oo)

    #hint_size = n*d #we truncate only t_A as opposed to [LNP22]

    #size_cut_tA   = n*d*(lgq-D) + (ell+lamb+1)*d*lgq+ m1*d*(log(G_entropy_const*sigma1, 2)) + m2*d*(log(G_entropy_const*sigma2_trunc,2))+ve*d*(log(G_entropy_const*sigmae,2))
    #size_cut_tA   = size_cut_tA + challnge_size + hint_size
    #size_BG_cut_tA   = n*d*(lgq-D) + (ell+lamb+1)*d*lgq+ m1*d*(log(G_entropy_const*sigma1, 2)) + (m2-n)*d*(log(G_entropy_const*sigma2_trunc,2))+ve*d*(log(G_entropy_const*sigmae,2))
    #size_BG_cut_tA   = size_BG_cut_tA + challnge_size + hint_size
    return size_plain/(8.*Kilo), size_cut/(8.*Kilo)


#####################
# Falcon parameters
# Adapted from 1st parameter set of the Falcon NIST submission
# Available at: https://falcon-sign.info/falcon.pdf
#####################
print(' ------------------ Commitment area  ------------------ ' )

d   = 128
sec = 128
n_F = 256
q_F = 1125899906842829

#####################
# Module-LWE security for c*u1+u2
# (replace the stat. arguments of LHL to computational LWE)
# h' is of the same dimension as Falcon's h
# and we keep the same q
#####################
tau_x = 2**22 #ell_oo norm of x1, x2
c_lwe_hardness = MLWE_optimize_attack(q_F, n_F, n_F, tau_x, cost_attack=LWE_primal_cost, cost_svp=svp_classical, verbose=False)
print('LWE hardness for c:', c_lwe_hardness)

print('Query size:', (256 * log(q_F, 2))/(8 * 1024))