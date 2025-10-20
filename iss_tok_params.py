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
from estimator_ntru import combined_attack_prob

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
print(' ------------------ Falcons  signature area  ------------------ ' )

d   = 128
sec = 128
n_F = 512
q_F = 7213
# This is not the modulus chosen by Falcon. We chose the smallest modulus that satisfies:
# (1) q_F mod 8 = 5 -> x^128+1 splits into two nice factors (see Lemma 2.5 in LNP)
# (2) q_F > 2^{128/lmbda}
# the smallness of q_F makes beta_inf small, hence milder lower bound on q_ZK later


epsinv = sqrt(sec*2**64)# R\'enyi divergence stuff, see above (2.13) in Falcon doc
#epsinv = 2**sec # stat. distance
sigma_F = (1/pi)*sqrt(log(4*n_F*(1+epsinv))/2)*1.17*sqrt(q_F)
print('sigma_F:', sigma_F)
# (2.13) in Falcon doc
beta = sigma_F*1.1*sqrt(2*n_F)
# (2.14) in Falcon doc, l2-norm of the signature
beta_inf = ceil(sigma_F*4.15)
# l2-inf norm bound of the sig, see [Le. 2.2, LNP22] with md=1.
# we chose that as large as possible to make the lower bound on encryption q_PKE to be as close to our q_PKE as possible
# Probability that a signature has the appropriate norm
# Maple code: t := 4.15: n:=512: (1-t*exp((-t^2 + 1)/2))^(n);
# Maple answer is 0.528

print('beta_inf:', beta_inf)
print('beta:', beta)

#Falcon's hardness (use the NTRU hardness estimator)
sigma_fg = 1.17*sqrt(q_F/(2*n_F)) # see Alg.5 in https://falcon-sign.info/falcon.pdf
print('sigma_fg:', sigma_fg)
beta_ntru_skr = combined_attack_prob(q_F, n_F, sigma_fg*sigma_fg, "circulant", 8, "SKR", verbose=False)
print(beta_ntru_skr[0], beta_ntru_skr[1])
print('Falcons secret key security as NTRU:', svp_classical(beta_ntru_skr[1]))

signature_key_primal = MLWE_optimize_attack(q_F, n_F, n_F, sigma_fg, cost_attack=LWE_primal_cost, cost_svp=svp_classical, verbose=False)
print('Falcons secret key security as LWE:', signature_key_primal)

# commented out because it takes quite some time and the attack is worse than primal
#forgery_primal = SIS_optimize_attack(q_F, 2*n_F, n_F, beta_inf, cost_attack=SIS_linf_cost, cost_svp=svp_classical, verbose=False)
#print('Falcons forgery as SIS in ell_infinity:', forgery_primal)

# returns 0 bits of security when beta>q_F
forgery_primal_ell2 = SIS_optimize_attack(q_F, 2*n_F, n_F, beta, cost_attack=SIS_l2_cost, cost_svp=svp_classical, verbose=False)
print('Falcons forgery as SIS in ell_2:', forgery_primal_ell2)

#####################
# ZK proof    (see LNP22, Fig. 10)
#####################

print(' ------------------ ZK proof  area  ------------------ ' )

# prove b, r are binary
ve = 2 # number of norm equations as proxy for proving (b * (b - 1) = 0 and r * (r - 1) = 0))
vd = 0 # number of inf-norm equations

norm_szk = q_F
# norm of the norm vector s^(e)
alpha_e = q_F
print('alpha_e:', alpha_e)

gamma_e = 3
t       = 1.64
kappa   = 2
Be = 2*sqrt(256/26)*t*gamma_e*sqrt(337)*(alpha_e) # as per Thm. 5.3
ce = d*(round(7*n_F/d)+1) # Fig 10
lb1 = floor(Be*41*ce)+1                      # Thm. 5.3
lb2 = floor(Be**2+sqrt(ve*d)*Be)+1           # Thm. 5.3
lb3 = floor(Be**2+2)  # Thm. 5.3
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

q2 = 2317417901 #chosen as next prime to ceil(lb/q) such that it is congruent to 5 mod 8
# Maple code: q3 := 124781: d:=128: l := numelems((Factors(x^d + 1) mod q3)[2]): isprime(q1), l;
# Maple answer is "true, 2"
q_ZK = q_F*q2
assert(q_ZK>=lb)
print('q3:', q2)
print('log(q_ZK):', log(q_ZK,2), 'q_ZK:', q_ZK)

l = 2 # number of factors of (x^128+1) mod q1

#summary of all the params
#  We are somewhat free to choose n and m2
#  For MLWE to make sense we require m2>n+{0,1,2}*(256/d+1)+lamb/2+ell,

gamma_candidate = 2**24
D_candidate = floor(log(gamma_candidate / (kappa*d), 2) + 1) # as per LNP22 (paragraph Dilithium compression)
print('D:', D_candidate)
iss_tok = {'d': 128, #ring dimension
              'n': 11, # determines SIS hardness
              'm1': 4, # m1
              'm2': 34, # m2 - n+{0,1,2}*(256/d+1)+lamb/2+ell determines LWE hardness
              'q': q_ZK, # ZK modulus
              'lambda': 10,
             'nu': 1,  # taken from LNP22
             'eta': 59, # taken from LNP22
             'gamma1': 10, # influences the repetition factor
             'gamma2': 1.5,  # influences the repetition factor
             'gammae': 5,  # influences the repetition factor
             'alphae':(alpha_e)**2, # norm of the norm vector s^(e) = [y1| y2 | x1 | x2 | s | e1 | e2 | x)
             'norm_s1': norm_szk**2, # norm of the zk secret
             'ell': 0,
             've': 3,  # number of exact norm proofs
             'vd': 0,  # number of approxumate norm proofs
             'D': D_candidate, # cutting D low-order bits from t_A
             'gamma': gamma_candidate # cutting log_gamma bits from w
             }
scalar = 1 #scalar is the number from the set {0,1,2} (see formula above) that depends on the number ve and vd
assert(iss_tok['m2']>iss_tok['n']+scalar*(256/d+1)+iss_tok['lambda']/2+iss_tok['ell']+2)

# sizes
proof = proof_size(iss_tok)
print('proof size:', proof)
#print('overall BG-optimized:', (ct_size+proof[1]))
#print('overall BG-optimized+CUT:', (ct_size+proof[3]))

# number of repetitions (rejection sampling)
# See [LNP22] Section 6.2  (12 is erroneous and should be 14)
reps = 2*exp(14/iss_tok['gamma1'] + 1/(2*iss_tok['gamma1']**2) + 1/(2*iss_tok['gamma2']**2) + 1/(2*iss_tok['gammae']**2))
print('avg nbr of repeats:', reps)

# security
sis_sec = SIS_security(iss_tok)
print('SIS hardness:', sis_sec)
lwe_sec = LWE_security(iss_tok)
print('LWE hardness primal:', lwe_sec)
#lwe_sec_dual = LWE_security(iss_tok, attack=LWE_dual_cost)
#print('LWE hardness dual:', lwe_sec_dual)
