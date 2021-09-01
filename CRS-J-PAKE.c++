#include <bits/stdc++.h>
#include <cmath>
#include <gmp.h> 
#include <gmpxx.h>
#include <iostream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <bitset>
#include <type_traits>
#include <string>
#include <openssl/sha.h>

using namespace std;

// Randomness functions 

unsigned long  getRandomInt(){ 
    std::random_device rd;
    return rd();
}

mpz_class random_number_range (mpz_class n) {
    gmp_randclass X(gmp_randinit_default) ; 
    X.seed(getRandomInt());
    mpz_class z = X.get_z_range(n);   
    return z ; 
}

//hashing algorithm 

string sha256(const string str)
{
    unsigned char hash[SHA256_DIGEST_LENGTH];
    SHA256_CTX sha256;
    SHA256_Init(&sha256);
    SHA256_Update(&sha256, str.c_str(), str.size());
    SHA256_Final(hash, &sha256);
    stringstream ss;
    for(int i = 0; i < SHA256_DIGEST_LENGTH; i++)
    {
        ss << hex << setw(2) << setfill('0') << (int)hash[i];
    }
    return ss.str();
}

//Ellyptic curves implementation 

struct Point 
{
    mpz_class x; 
    mpz_class y; 
};

struct Elliptic_Curve
{
    mpz_class p;
    mpz_class a; 
    mpz_class b;
    struct Point g; 
    mpz_class n;
    mpz_class h; 
};

Elliptic_Curve initialize_curve() {
    Elliptic_Curve ECC; 
    ECC.p = ("115792089237316195423570985008687907853269984665640564039457584007908834671663");
    ECC.a = ("0");
    ECC.b = ("7");
    ECC.n = ("115792089237316195423570985008687907852837564279074904382605163141518161494337");
    ECC.h = ("1");
    ECC.g.x = ("55066263022277343669578718895168534326250603453777594175500187360389116729240");
    ECC.g.y = ("32670510020758816978083085130507043184471273380659243275938904335757337482424");
    return ECC; 
}

bool is_on_curve(Point a, Elliptic_Curve curve) {
    if (a.x == 0 && a.y == 0) {
        return true; 
    }
    mpz_class x = a.x; 
    mpz_class y = a.y; 

    if (((y * y - x * x * x  - curve.a * x - curve.b) % curve.p) == 0) {
        return true; 
    }

    else {
        return false; 
    }
}


mpz_class inverse_mod(mpz_class k, mpz_class p) {
    mpz_class temp ; 
    if (k == 0) {
        throw "Division by zero condition!";
    }

    if (k < 0) {
        return p - inverse_mod(-k, p); 
    }

    mpz_class s = 0 ; 
    mpz_class old_s = 1; 
    mpz_class z = 1;  
    mpz_class old_z = 0; 
    mpz_class r = p;
    mpz_class old_r = k; 

    while (r!= 0) {
        mpz_class quotient; 
        mpz_fdiv_q(quotient.get_mpz_t(), old_r.get_mpz_t(), r.get_mpz_t()); 
        temp = old_r - quotient * r; 
        old_r = r; 
        r = temp; 
        temp = old_s - quotient * s; 
        old_s = s; 
        s = temp; 
        temp = old_z - quotient * z; 
        old_z = z; 
        z = temp;  
    } 

    mpz_class gcd = old_r; 
    mpz_class x = old_s; 
    mpz_class y = old_z; 

    assert (gcd == 1);
    return x % p ;  
}


Point point_neg(Point a, Elliptic_Curve curve) {
    assert (is_on_curve(a,curve)); 

    if(a.x == 0 && a.y == 0) {
        return a; 
    }

    Point result; 
    result.x = a.x; 
    result.y = -a.y % curve.p;

    assert (is_on_curve(result, curve)); 

    return result;  
}

Point point_add(Point point1, Point point2, Elliptic_Curve curve) {
    assert (is_on_curve(point1, curve));
    assert (is_on_curve(point2, curve)); 

    if (point1.x == 0 && point1.y == 0) {
        return point2; 
    }
    if (point2.x == 0 && point2.y == 0) {
        return point1; 
    }
    mpz_class x1 = point1.x; 
    mpz_class y1 = point1.y; 

    mpz_class x2 = point2.x; 
    mpz_class y2 = point2.y; 

    if (x1 == x2 && y1 != y2) {
        Point point0; 
        point0.x = 0; 
        point0.y = 0; 
        return point0; 
    }

    mpz_class m; 

    if (x1 == x2) {
        m = (3 * x1 * x1 + curve.a) * inverse_mod(2*y1, curve.p);
    }
    else {
        m = (y1 - y2) * inverse_mod(x1 -x2, curve.p); 
    }

    mpz_class x3 = m*m - x1 -x2 ; 
    mpz_class y3 = y1 + m * (x3 - x1); 
    Point result; 
    result.x = x3 % curve.p; 
    result.y = -y3 % curve.p; 

    assert (is_on_curve(result, curve)); 

    return result; 
}

void CSWAP(mpz_class b, mpz_class *x, mpz_class *y) {
    mpz_class m = 0 - b; 
    mpz_class z = (*x ^ *y);
    mpz_class d = m & z; 
    mpz_class xp = *x ^ d; 
    mpz_class yp = *y ^ d; 
    *x = xp; 
    *y = yp;  
}

mpz_class reduction(mpz_class a, Elliptic_Curve curve) {
    mpz_class u; 
    u = ("463168356949264781694283940034751631414809620208824894785240019497232391142146"); 
    mpz_class t = (a * u);  
    t = t >> 514;
    t = a - (t * curve.n); 
    mpz_class b = (((curve.n+((~t)+1)) >> 514) & 1);
    u = t; 
    t = t - curve.n; 
    CSWAP(b, &t, &u);
    return u; 
}


Point scalar_mult(mpz_class k, Point point, Elliptic_Curve curve) {
    assert (is_on_curve(point, curve)); 

    Point R0; 
    R0.x = 0; 
    R0.y = 0;

    if ((k % curve.n) == 0 || (point.x == 0 && point.y == 0)) { 
        return R0; 
    }

    if (k<0) {
        return scalar_mult(-k, point_neg(point,curve), curve); 
    }

    k = reduction(k, curve); 
    int size = mpz_sizeinbase(k.get_mpz_t(), 256);   // 256-bit = 32-byte
    unsigned char scalarbyte[size] = { 0 };          // initialize scalar byte array 
    // export scalar k to 32-byte array
    mpz_export(scalarbyte, NULL, -1, 1, 1, 0, k.get_mpz_t());
    
    Point R1; 
    R1.x = point.x;
    R1.y = point.y;  

    // Montgomery ladder 
    for (int i = 255; i >= 0; i--) {
        int sbit = (scalarbyte[i>>3] >> (i & 7)) & 1; // obtain the current scalar bit 
        int sbit1 = (scalarbyte[(i+1)>>3] >> ((i+1) & 7)) & 1; 
        int b = sbit ^ sbit1;
        CSWAP(b, &R0.x, &R1.x);  
        CSWAP(b, &R0.y, &R1.y);
        R1 = point_add(R0, R1, curve); 
        R0 = point_add(R0, R0, curve);
    }
    
    int sbit = (scalarbyte[0>>3] >> (0 & 7)) & 1; 
    CSWAP(sbit, &R0.x, &R1.x);  
    CSWAP(sbit, &R0.y, &R1.y);
    assert (is_on_curve(R0,curve)); 
    return R0; 
}

mpz_class hash_to_int(string s){
    std::stringstream ss;
    mpz_class x; 
    s = sha256(s); 
    ss << std::hex << s;
    ss >> x; 
    return x; 
}

Point generation_CRS(Elliptic_Curve curve) {
    mpz_class r = random_number_range((curve.n-1));
    Point u = scalar_mult(r, curve.g, curve); 
    return u; 
}

string point_to_string(Point point) {
    std::stringstream ss;
    string result1;
    string result2;  
    ss << point.x; 
    ss >> result1;
    std::stringstream().swap(ss); 
    ss << point. y; 
    ss >> result2; 
    return result1+result2; 
}

mpz_class point_to_mpz(Point point) {
    string str = point_to_string(point); 
    std::stringstream ss;
    ss << str; 
    mpz_class m ; 
    ss >> m; 
    return m;  
} 

string stm(mpz_class x) {
    string str; 
    std::stringstream ss;
    ss << x;
    ss >> str;
    return str; 
}

bool check_two_points(Point a, Point b) {
    if (a.x == b.x) {
        return true; 
    }
    else {
        return false; 
    }
}

void generate_nizk_proof(Point G, mpz_class x, Point XS, Point ID, string U, Elliptic_Curve curve, Point *r1, mpz_class *r2) {
    mpz_class v = random_number_range(curve.n); 
    Point vS = scalar_mult(v, G, curve); 
    string chal = point_to_string(G) + point_to_string(XS) + point_to_string(ID) + U + point_to_string(vS); 
    mpz_class c = hash_to_int(chal); 
    mpz_class pi = reduction((v-c*x), curve); 
    *r1 = vS; 
    *r2 = pi; 
}

void generate_nizk_proof_2(Point G, mpz_class x, Point XS, string L, Elliptic_Curve curve, Point *r1, mpz_class *r2) {
    mpz_class v = random_number_range(curve.n); 
    Point vS = scalar_mult(v, G, curve); 
    string chal = point_to_string(G) + point_to_string(XS) + L + point_to_string(vS); 
    mpz_class c = hash_to_int(chal); 
    mpz_class pi = reduction((v-c*x), curve); 
    *r1 = vS; 
    *r2 = pi;
}

void check_nizk_proof(Point G, mpz_class pi, Point X, Point ID, string U, Point Vs, Elliptic_Curve curve) {
    string chal = point_to_string(G) + point_to_string(X) + point_to_string(ID) + U + point_to_string(Vs); 
    mpz_class c = hash_to_int(chal); 
    Point check = point_add(scalar_mult(pi,G,curve), scalar_mult(c, X,curve), curve); 
    if (check_two_points(check, Vs)) {
        cout << ("It has been proved") << endl; 
    }
    else {
        cout << ("not proved") << endl;
    }
}

void check_nizk_proof_2(Point G, mpz_class pi, Point X, string  L, Point Vs, Elliptic_Curve curve) {
    string chal = point_to_string(G) + point_to_string(X) + L + point_to_string(Vs); 
    mpz_class c = hash_to_int(chal); 
    Point check = point_add(scalar_mult(pi,G,curve), scalar_mult(c, X,curve), curve); 
    if (check_two_points(check, Vs)) {
        cout << ("It has been proved") << endl; 
    }
    else {
        cout << ("not proved") << endl;   
    }
}

int main() {

    // Initialization

    Elliptic_Curve curve = initialize_curve(); 
    double duration;
    std::clock_t start;
    start = std::clock();
    string secret = "IssamJomaatest"; 
    mpz_class s = hash_to_int(secret); 

    Point G = curve.g; 
    Point U = generation_CRS(curve);
    string strU = point_to_string(U); 

    //Round1

    mpz_class x1 = random_number_range(curve.n); 
    mpz_class x2 = random_number_range(curve.n); 

    Point X1 = scalar_mult(x1, curve.g, curve);
    Point X2 = scalar_mult(x2, curve.g, curve); 

    Point A = point_add(X1, X2, curve); 
    A = scalar_mult(x1*s, A, curve);

    Point B = point_add(X1, X2, curve); 
    B = scalar_mult(x2*s, A, curve); 

    Point Vs; 
    mpz_class pi_1; 
    generate_nizk_proof(G, x1, X1, A, strU, curve, &Vs, &pi_1);
    check_nizk_proof(G, pi_1, X1, A, strU, Vs, curve); 
    
    Point Vg; 
    mpz_class pi_2; 
    generate_nizk_proof(G, x2, X2, B, strU, curve, &Vg, &pi_2); 
    check_nizk_proof(G, pi_2, X2, B, strU, Vg, curve); 

    //Round2

    Point beta = point_add(U, X1, curve); 
    beta = scalar_mult(x2*s, beta, curve); 

    Point alpha = point_add(U, X2, curve); 
    alpha = scalar_mult(x1*s, alpha, curve); 

    Point UX1 = point_add(U, X1, curve); 
    Point UX2 = point_add(U, X2, curve); 

    string la = point_to_string(A) + point_to_string(B) + point_to_string(X1) + point_to_string(X2) + strU; 
    string lb = point_to_string(B) + point_to_string(A) + point_to_string(X2) + point_to_string(X1) + strU; 

    Point Valpha;
    mpz_class pi_alpha; 
    generate_nizk_proof_2(UX2, x1*s, alpha, la, curve, &Valpha, &pi_alpha); 
    check_nizk_proof_2(UX2, pi_alpha, alpha, la, Valpha, curve); 

    Point Vbeta; 
    mpz_class pi_beta; 
    generate_nizk_proof_2(UX1, x2*s, beta, lb, curve, &Vbeta, &pi_beta); 
    check_nizk_proof_2(UX1, pi_beta, beta, lb, Vbeta, curve); 

    //Final Round

    Point ka = scalar_mult(-x1*s, X2, curve);
    ka = point_add(beta, ka, curve); 
    ka = scalar_mult(x1, ka, curve); 

    Point kb = scalar_mult(-x2*s, X1, curve); 
    kb = point_add(alpha, kb, curve); 
    kb = scalar_mult(x2, kb, curve); 

    mpz_class ska = hash_to_int(stm(ka.x));
    mpz_class skb = hash_to_int(stm(kb.x));
    
    duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
    cout << duration << endl; 
    cout << "This is the result of the Client A key : " << ska << endl; 
    cout << "This is the result of the Client B key : " << skb << endl;

    if (ska == skb) {
        cout << "This test was a success." << endl; 
    }
}