from sympy import *
from sympy.polys.galoistools import *
from sympy.polys.factortools import *
x = S('x')
 # for mod 2 case

def trace(x, d):
    s = 0
    for i in range(d):
        s += x**(2**i)
    return expand(s)
####EQF test###

def EQ(v1,i,L,d,mv,count = 1):
    
    if mv != 2 :
        if degree(v1) != 1 :
            #print("check",degree(v1))
            if True : # must be odd number
            #d = 2
            #i = 1
            #mv = 11
                #print("deg mv", mv)
                md =(mv**d - 1)/2
                qq = (x + i)**int(md)
                print(qq)
                bound1 = degree(v1)/d
                secfac = v1
                
            #for count in range(1,int(bound1)):
                #print('\n Loop \n')
                tmp1 = gcd(secfac,qq - 1,modulus = mv  , symmetric = True)
                
                    
                while tmp1 == 1 :
                    i = i + 1
                    qq = (x + i)**md
                    tmp1 = gcd(secfac,qq - 1 ,modulus = mv , symmetric = True)
                    print(tmp1)
                       
                #print("tmp1 = ", tmp1)
                if degree(tmp1)>=2 and not degree(tmp1)==d:
                    EQ(tmp1,i+1,L,d,mv)
                else :
                
                    L.append(tmp1)
                secfac = quo(secfac,tmp1,modulus = mv , symmetric = True)
                
                #print("second factor :",secfac)
                
                if degree(secfac)>=2 and not  degree(secfac)==d :
                    EQ(secfac,i+1,L,d,mv)
                else :
                    L.append(secfac)
                
                i =  i + 1
                qq = (x + i)**md
            #print('L is ',L)
            count= count + 1
            return L
    else:
        
        d1 = d
        if count == 1:
            #print('check1')
            r = Poly(v1).degree_list()[0]/d
            rd = r*d - 1
        else:
            #print(check)
            r = Poly(v1).degree_list()[0]
            rd = r - 1
        i = 0
        qqq = x**int(rd) + 1
        #print('q ', qqq)
        secfac = v1
        #tr = trace(qqq,d)
        tr = trace(x,d1)
        expr1 = tr.subs({x : qqq})

        poly = expr1.as_poly(x).set_modulus(2)
        #print('polytrace ', poly)
        tr = poly.as_expr(x)
        #print('trace ',tr)
        tmp1 = gcd(secfac,tr,modulus = mv , symmetric = True)
        c = 0
        l = int(rd) 
        while (tmp1 == 1 or tmp1 == 0) and l >= 1:
             
            i = i + 1
            #print(i,l)
            if i == 1 :
                qqq = x**int(rd) + x**i
            
            else :
                l = l - 1
                if l!=1:
                    #print(x**l,'   ',qqq)
                    qqq = x**int(rd) + x**(l) + 1
            #print('qqq ',qqq)   
            #qqq  = qqq + x**i
            tr = trace(x,d1)
            expr1 = tr.subs({x : qqq})

            poly = expr1.as_poly(x).set_modulus(2)
            #print('polytrace ', poly)
            tr = poly.as_expr(x)
            #print('trace ',tr)
            tmp1 = gcd(secfac,tr,modulus = mv , symmetric = True)
            #print('tmp1 ',tmp1)
        #print('tmp1',tmp1)
        if tmp1 == 1 :
        
            
            qqq = x**int(rd)  
            tr = trace(x,d1)
            expr1 = tr.subs({x : qqq})

            poly = expr1.as_poly(x).set_modulus(2)
            #print('polytrace ', poly)
            tr = poly.as_expr(x)
            #print('trace ',tr)
            tmp1 = gcd(secfac,tr,modulus = mv , symmetric = True)
            #print('tmp1 ',tmp1)
        
        
        
        
        
        
        
        if not degree(tmp1) == d1 :
            count = count + 1
            EQ(tmp1,1,L,d,mv)
        else :    
            L.append(tmp1)
        secfac = quo(secfac,tmp1,modulus = mv , symmetric = True)
        #print('secfac ',secfac)
        if secfac != 1 :
            if  not  degree(secfac)==d1 :
                    count = count + 1
                    EQ(secfac,1,L,d,mv)
            else :
                    L.append(secfac)
        
            
    #print('L is ',L)        
    return L     
            
            
