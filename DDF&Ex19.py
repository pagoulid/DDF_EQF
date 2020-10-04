from sympy import *
from sympy.polys.galoistools import *
from sympy.polys.factortools import *
from EQF import EQ
x = S('x')








def DDFactor(p,m):
    pi = []

    r = 0
    old = 0


    while True :
        r = r+1
        exp = m**r
        q = x**int(exp) - x
    
        if r == 1:
            v = gcd(p,q,modulus = m , symmetric = True)
            
        
        elif r == 2 :
            old = quo(p,pi[0], modulus = m , symmetric = True)
            v=gcd(old,q , modulus = m , symmetric = True)
            
        else:
            
            old = quo(old,pi[r-2] , modulus = m , symmetric = True) # p/(p1...pn) 
        
            v=gcd(old,q, modulus = m , symmetric = True)
            
        if old == 1 :
            break
        else :
            
            pi.append(v)
        
        #print("Factorization of {} mod {} is : {}".format(p,m,pi))
    print(pi)
    if 1 in pi:
        pi=list(filter(lambda a: a != 1, pi))
        
        
    factors = ''
    for p in pi :
        factors = factors + '*' + '( '+str(p) + ' )'
    return factors[1:],pi






x = S('x')
#p = x**15 + 1     # x**12 + x**9 + x**6 + x**3 + 1 m = 2 d = 4
#p = x**5 - 12*x**3 + 7*x**2 + 35*x - 35 # x**2 - 5  m = 11 d = 1
#p = x**4 + 5*x**3 + 6*x**2 + 5*x + 1 m = 5 d = 2


#p = x**10 + x**5 + 1 # x**8 + x**7 + x**5 + x**4 + x**3 + x + 1 d = 4 m = 2
#p = x**5 + x + 1
#p = x**4 - 8*x - 9
#p = x**5 + 1

def prime(p):

    i = 2
    nxt = 2
    sprime=gf_sqf_p( Poly(p).all_coeffs(), nxt, ZZ )
    print('prime : {} , is suitable : {}'.format(nxt,sprime))
    
    while sprime == False :
        #i = i+1
        nxt = nextprime(nxt,ith=i)
        sprime=gf_sqf_p( Poly(p).all_coeffs(), nxt, ZZ )
        print('prime : {} , is suitable : {}'.format(nxt,sprime))

    print('m =',simplify(nxt))
    m = nxt
    return m


  
  
  
  



p = x**7 + 1
m= prime(p)



new,h=DDFactor(Poly(p, x, domain=GF(m, symmetric=True)).as_expr(),m)
print("Factorization of {} mod {} is : {}".format(p,m,new))
d = 3
v1 = x**6 + x**5 + x**4 + x**3 + x**2 + x + 1
d = 3
L = EQ(v1,0,[],d,m)
q=factor(p, modulus = m,symmetric = True).as_ordered_factors()

#cont = h._hashable_content()
print(L)
print(q)
print('\n ##################################### \n')
print('\n')
print("Factorization of {} mod {} is : {}".format(p,m,new))
print('\n Final : (',h[0],')*','(',L[0],')*(',L[1],')')
print('evaluate :',factor(p,modulus = m , symmetric = True))
print('\n')
print('\n ################# \n')




p = x**21 + 1
m = prime(p)

new,h=DDFactor(Poly(p, x, domain=GF(m, symmetric=True)).as_expr(),m)

v1 = x**6 + x**5 + x**4 + x**3 + x**2 + x + 1
# x**12 + x**11 + x**9 + x**8 + x**6 + x**4 + x**3 + x d = 6
d = 3

L1 = EQ(v1,0,[],d,m) # solution



v2 = x**12 + x**11 + x**9 + x**8 + x**6 + x**4 + x**3 + x + 1
d = 6

L2 = EQ(v2,0,[],d,m) # solution
print('\n ################# \n')
print('\n')
print("Factorization of {} mod {} is : {}".format(p,m,new))
print('\n Final : (',h[0],')*(',h[1],')*(',L1[0],')*(',L1[1],')*(',L2[0],')*(',L2[1],')')
print('evaluate :',factor(p,modulus = m , symmetric = True))
print('\n')
print('\n ################# \n')


p = x**22 + x**21 + x**20 + x**19 + x**18 + x**17 + x**16 + x**7 + x**6 + x**5 + x**4 + x**3 + x**2 + x
m = prime(p)
new,h=DDFactor(Poly(p, x, domain=GF(m, symmetric=True)).as_expr(),m)

v3 = x**6 + x**5 + x**4 + x**3 + x**2 + x + 1
d = 3

L3 = EQ(v3,0,[],d,m) # solution

v4 =x**12 + x**9 + x**6 + x**3 + 1
d = 4
L4 = EQ(v4,0,[],d,m) # solution

print('\n ################# \n')
print('\n')
print("Factorization of {} mod {} is : {}".format(p,m,new))
print('\n Final : (',h[0],')*(',h[1],')*(',L3[0],')*(',L3[1],')*(',L4[0],')*(',L4[1],')*(',L4[2],')')
print('evaluate :',factor(p,modulus = m , symmetric = True))
print('\n')
print('\n ################# \n')
