from secp256k1 import curve,scalar_mult, point_add,point_neg
import random
import hashlib
import sys
import time

def generation_CRS(): #Generation of the CRS string I doubt that this is a proper generation method. 
    r = random.randint(0, curve.n-1)
    u = scalar_mult(r, curve.g)
    del r 
    return u 

def generate_nizk_proof(G, x, XS, ID, U):
    v = random.randint(0, curve.n-1)
    vS = scalar_mult(v, G)
    chal = str(G)+str(XS)+str(ID)+str(U)+str(vS)
    c = int(hashlib.sha256(chal.encode()).hexdigest(),16)
    pi = (v-c*x) % (curve.n)
    return vS, pi 
    
def generate_nizk_proof_2(G, x, XS, L):
    v = random.randint(0, curve.n-1)
    vS = scalar_mult(v, G)
    chal = str(G)+str(XS)+str(L)+str(vS)
    c = int(hashlib.sha256(chal.encode()).hexdigest(),16)
    pi = (v-c*x) % (curve.n)
    return vS, pi 

def check_nizk_proof(G, pi, X, ID, U, Vs):
    chal = str(G)+str(X)+str(ID)+str(U)+str(Vs)
    c = int(hashlib.sha256(chal.encode()).hexdigest(),16)
    check = point_add(scalar_mult(pi,G), scalar_mult(c, X))
    if (check == Vs) : 
        print ("It has been proved")
    else : 
        print("not proved")
        quit()

def check_nizk_proof_2(G, pi, X, L, Vs):
    chal = str(G)+str(X)+str(L)+str(Vs)
    c = int(hashlib.sha256(chal.encode()).hexdigest(),16)
    check = point_add(scalar_mult(pi,G), scalar_mult(c, X))
    if (check == Vs) : 
        print ("It has been proved")
    else : 
        print('not proved')
        quit()

# Sorry for the variable names I had no inspiration. 
start = time.time()

secret="IssamJomaatest"
s=int(hashlib.sha256(secret.encode()).hexdigest(),16)


########################## Part 1 Generating variables.###############################
   
G = curve.g
U = generation_CRS() 


#Client A generate x1
x1 = random.randint(0, curve.n-1)
#Server B generate x2
x2 = random.randint(0, curve.n-1)

#Client A generate X1
X1 = scalar_mult(x1, curve.g)
#Server B generate X2
X2 = scalar_mult(x2, curve.g)

#A = (X1 X2) * (x1*s) not sure about the generation of this variable there is no indication on the J-pake-CRS paper so 
# I improvised based on the J-pake example I had. 

A = point_add(X1,X2)
A = scalar_mult(x1*s, A)

#B = (X1+X2) * (x2*s) same thing as A
B = point_add(X1,X2)
B = scalar_mult(x2*s, A)


############################### Part 2, two first NIZK proofs ############################


if (X1 == 1 ): 
        print("the protocol is aborted there is an error")
        quit()
if (X1 == 1 ): 
        print("the protocol is aborted there is an error")
        quit()


# Client A generate the first proof and "Theoretically send Vs(The random value created during the generation of the NIZK), 
# pi_1(equals to v-c*x), and A which is the ID of the user"
Vs,pi_1 = generate_nizk_proof(G,x1,X1,A,U)

# Server B once he received the variables will check the NIZK proof, 
check_nizk_proof(G,pi_1,X1,A,U,Vs)

#Now it's the turn of Server B to generate a proof 
Vg,pi_2 = generate_nizk_proof(G,x2,X2,B,U)

#server A will now verify the proof 
check_nizk_proof(G,pi_2, X2, B, U, Vg)


################################################ Part 3: another round of variable generation ########################
# beta = (U X1) * (x2 pw)
beta = point_add(U, X1)
beta = scalar_mult(x2*s, beta)

# alpha = (U X2) * (x1 pw)
alpha = point_add(U, X2)
alpha = scalar_mult(x1*s, alpha)

#gezneration of UX1
UX1 = point_add(U, X1)
#generation of UX2
UX2 = point_add(U, X2)

#Generation if the labels 
la = str(A)+str(B)+str(X1)+str(X2)+str(U)
lb = str(B)+str(A)+str(X2)+str(X1)+str(U)


################################################ Part 4: two last NIZK proofs #####################################
#Generation of the third NIZK by Client A to prove x1s 
Valpha, pi_alpha = generate_nizk_proof_2(UX2, x1*s, alpha, la)
#Server B Check the 2nd Nizk proof 
check_nizk_proof_2(UX2, pi_alpha, alpha, la, Valpha)

#Generation of the forth NIZK by the Server B to prove x2s 
Vbeta, pi_beta = generate_nizk_proof_2(UX1, x2*s, beta, lb)
#Client A check the forth NIZK 
check_nizk_proof_2(UX1, pi_beta, beta, lb, Vbeta)



############################################## Last variable generation and fiding the  session key #####################

# ka = (α X1*(−x2 pw)) * x2
ka = scalar_mult(-x1*s, X2)
ka = point_add(beta, ka)
ka = scalar_mult(x1, ka)

# kb = (β X2*(−x1 pw)) * x1
kb = scalar_mult(-x2*s, X1) 
kb = point_add(alpha, kb)
kb = scalar_mult(x2, kb)

# last step the results are converted into strings and hashed using SHA256. 
ska=int(hashlib.sha256(str(ka[0]).encode()).hexdigest(),16)
skb=int(hashlib.sha256(str(kb[0]).encode()).hexdigest(),16)

end = time.time()

print("this is the result of the Client A key : ", ska, "\n")
print("this is the result of the Server B key : ", skb, "\n")

if (ska == skb) : 
    print("The Test was a success")

print(end - start)

print(scalar_mult(8,curve.g))