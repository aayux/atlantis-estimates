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
q_F = 7213

#####################
# Module-LWE security for c*u1+u2
# (replace the stat. arguments of LHL to computational LWE)
# h' is of the same dimension as Falcon's h
# and we keep the same q
#####################
tau_x = 10 #ell_oo norm of x1, x2
c_lwe_hardness = MLWE_optimize_attack(q_F, n_F, n_F, tau_x, cost_attack=LWE_primal_cost, cost_svp=svp_classical, verbose=False)
print('LWE hardness for c:', c_lwe_hardness)

#####################
# Encryption Parameters
#
# Use Module-LWE (Kyber) over Rq = Zq[x]/(x^128+1)
# we encrypt:
# -- y1 (the first part of Falcon's signature)
# -- x1, x2 (user's randomizers)
# We have 3 elements over Zq[x]/(x^512+1) to encrypt
#
# a1 \in Rq^{8 X 8}
# a2 = esp1 * a1 + eps2;   eps1, eps2 \in Rq^{12 X 8}
# ct1 = a1*s+e1;           s, e1 \in Rq^{8 X 1}
# ct2 = a2*s+e2+p*(s1|x1|x2);      e2, s1 \in Rq^{12 X 1}
#
######################

print(' ------------------ PKE  area  ------------------ ' )

hardness = 8 # set minimal to get sufficient hardness
n_Enc = d * hardness
tau = 3     # inf-norm bound on r_i, e_i, eps1 and eps2

# Decryption must hold for ri's, ei's that can be as large as
# guaranteed by the ZK proof, i.e.,
# ||ri,ei1, ei2|| <= tau*sqrt(2*n_Enc+2*n_F)
# Independently, the sk eps1, eps2 is generated by the
# challenger and satisfies ||eps1,eps_2||_oo <= tau.
# We must have || eps1*s+e2-eps2*e1 ||_oo < p/2
# same for the randomness of ct2
boundp= sqrt(2*n_Enc)*tau**2*sqrt(2*n_Enc+2*n_F) + tau*sqrt(2*n_Enc+2*n_F)
p = ceil(2*boundp)
print('p:', p)

# We must have || eps1*r1+e2_1-eps2*e1_1 + p(y1||x1||x2) ||_oo < q_enc/2
boundq = ceil(2 * (p/2. + p * tau_x))

# We set q_enc = q_F * q1, with four conditions:
# (0) q1 prime
# (1) lamb*log(q1, 2) >= 128
# (2) (x^d+1) has two factors mod q1
# (3) q_enc > boundq
# We take q1 = q_F
lamb = 10
lgq1 = 128./lamb
q1   = 7213 # 7213 % 8 = 5
# Maple code: q1 := 7187: d:=128: l := numelems((Factors(x^d + 1) mod q1)[2]): isprime(q1), l;
# Maple answer is true, 2
#
assert(log(q1)/log(2)>=lgq1)

q_PKE = q_F *  q1
print('lower boundq:', boundq)
print('q_PKE:', q_PKE)

assert(q_PKE>=boundq), 'encryption q_PKE is too small'


# Hardness of LWE instance from encryption
lwe_enc_primal = MLWE_optimize_attack(q_PKE, n_Enc, n_Enc+2*n_F, tau, cost_attack=LWE_primal_cost, cost_svp=svp_classical, verbose=False)
print('LWE encryption primal:', lwe_enc_primal)
#lwe_enc_dual = MLWE_optimize_attack(q_PKE, n_Enc, n_Enc+512, tau, cost_attack=LWE_dual_cost, cost_svp=svp_classical, verbose=False)
#print('LWE encryption dual:', lwe_enc_dual)

#####################
# ZK proof    (see LNP22, Fig. 10)
# zk_secret s_ZK = [u1 | u2 |  t | s | e1 | e2]
#####################

print(' ------------------ ZK proof  area  ------------------ ' )

ve = 2 # number of norm equations
vd = 0 # number of inf-norm equations
norm_x = sqrt(ve)*sqrt(d)

# norm of [u1 & u2 & t & s]
norm_szk = sqrt(3*(tau_x^2*n_F) + tau**2*(2*n_Enc+2*n_F))
#norm of the norm vector s^(e) = [u1| u2 | t | s | e1 | e2 | x)
alpha_e = sqrt(3*tau_x**2*n_F + tau**2*(2*n_Enc+2*n_F)  + norm_x**2)
print('alpha_e:', alpha_e)

gamma_e = 3
t       = 1.64
kappa   = 2
Be = 2*sqrt(256/26)*t*gamma_e*sqrt(337)*(alpha_e) # as per Thm. 5.3
ce = d*(round(2*n_F/d)+1+round((2*n_Enc+2*n_F)/d)+1) # Fig 10
lb1 = floor(Be*41*ce)+1                      # Thm. 5.3
lb2 = floor(Be**2+sqrt(ve*d)*Be)+1           # Thm. 5.3
beta1 = sqrt(2*n_F)*tau_x                    # ||(x1,x2)||
beta2 = tau*sqrt(n_Enc+n_Enc+12*128)         # ||(s,e1,e2)||
lb3 = floor(Be**2+2*(max(beta1, beta2))**2)  # Thm. 5.3
lb = max(lb1, lb2, lb3)


print(2*sqrt(256/26)*t*gamma_e*sqrt(337))
print('ce:', ce, 'log(Be):', log(Be,2))
print()
print('lb1:', log(lb1,2), 'lb2:', log(lb2,2), 'lb3:',log(lb3,2))

# q_ZK = q_PKE * q3  = q_F * q1 * q3
# We want
# (0) q3 prime
# (1) q3 >= q1 (q1=q_F should remain the smallest factor of q_ZK)
# (2) q_ZK>= lb
# (3) (x^d+1) has two factors mod q_ZK

q3 = 7213 #chosen as next prime to ceil(lb/q) such that it is congruent to 5 mod 8
# Maple code: q3 := 124781: d:=128: l := numelems((Factors(x^d + 1) mod q3)[2]): isprime(q1), l;
# Maple answer is "true, 2"
q_ZK = q_PKE*q3
assert(q_ZK>=lb)
assert(q3>=q1)
print('q3:', q3)
print('log(q_ZK):', log(q_ZK,2), 'q_ZK:', q_ZK)

l = 2 # number of factors of (x^128+1) mod q1
kappa_bound = (1./(2*sqrt(l)))*(q1)**(1./l)
# The sigma_{-1} row of Fig 3 can be used if kappa < kappa_bound
assert(kappa < kappa_bound)
print('kappa bound:', kappa_bound)

#summary of all the params
#  We are somewhat free to choose n and m2
#  For MLWE to make sense we require m2>n+{0,1,2}*(256/d+1)+lamb/2+ell,

gamma_candidate = 2**24
D_candidate = floor(log(gamma_candidate / (kappa*d), 2) + 1) # as per LNP22 (paragraph Dilithium compression)
print('D:', D_candidate)
cli_que = {'d': 128, #ring dimension
              'n': 13, # determines SIS hardness
              'm1': 16, # m1 = len([u1|u2|t]) + len(s|e1|e2) + ve = (3 * 2 + 2 * 2 + 4 + 2) = 16 (here s is the encryption secret)
              'm2': 36, # m2 - n+{0,1,2}*(256/d+1)+lamb/2+ell determines LWE hardness
              'q': q_ZK, # ZK modulus
              'lambda': lamb,
             'nu': 1,  # taken from LNP22
             'eta': 59, # taken from LNP22
             'gamma1': 10, # influences the repetition factor
             'gamma2': 1.5,  # influences the repetition factor
             'gammae': 5,  # influences the repetition factor
             'alphae':(alpha_e)**2, # norm of the norm vector s^(e)
             'norm_s1': norm_szk**2, # norm of the zk secret
             'ell': 0,
             've': 3,  # number of exact norm proofs
             'vd': 0,  # number of approxumate norm proofs
             'D': D_candidate, # cutting D low-order bits from t_A
             'gamma': gamma_candidate # cutting log_gamma bits from w
             }
scalar = 1 #scalar is the number from the set {0,1,2} (see formula above) that depends on the number ve and vd
assert(cli_que['m2']>cli_que['n']+scalar*(256/d+1)+cli_que['lambda']/2+cli_que['ell']+2)

# sizes
proof = proof_size(cli_que)
print('proof size:', proof)
ct_size = (n_Enc+2*n_F)*log(q_PKE,2)/(8.*Kilo)
print('ct size:', ct_size)
print('overall non optimized:', (ct_size+proof[0]))
print('overall with CUT:', (ct_size+proof[1]))
#print('overall BG-optimized:', (ct_size+proof[1]))
#print('overall BG-optimized+CUT:', (ct_size+proof[3]))

# number of repetitions (rejection sampling)
# See [LNP22] Section 6.2  (12 is erroneous and should be 14)
reps = 2*exp(14/cli_que['gamma1'] + 1/(2*cli_que['gamma1']**2) + 1/(2*cli_que['gamma2']**2) + 1/(2*cli_que['gammae']**2))
print('avg nbr of repeats:', reps)

# security
sis_sec = SIS_security(cli_que)
print('SIS hardness:', sis_sec)
lwe_sec = LWE_security(cli_que)
print('LWE hardness primal:', lwe_sec)
#lwe_sec_dual = LWE_security(cli_que, attack=LWE_dual_cost)
#print('LWE hardness dual:', lwe_sec_dual)

