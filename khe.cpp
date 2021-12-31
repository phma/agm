/******************************************************/
/*                                                    */
/* khe.cpp - compute khe function                     */
/*                                                    */
/******************************************************/

/* The khe function խ(z) is defined on the left half-plane as follows:
 * խ(z) is asymptotic to 4*exp(z)+1 as Re(z)-> -∞.
 * That is, lim[Re(z)-> -∞] (խ(z)-1)/exp(z) =4.
 * խ(z)*խ(z+πi)=խ(2z+πi)²
 * խ(z)+խ(z+πi)=2խ(2z)
 * խ(z) is continuous on the left half-plane.
 *
 * The function խ(x+yi), as a function of y where x is fixed, is called a loop
 * (it's periodic). To compute a loop, this program starts with 36 points on
 * a circle of radius 65 ulps centered on p in [?,2), where the circle centered
 * at p=1 is the loop for x=-ln(2^54/65) (-33.25556048034141) within machine
 * precision, assuming 8-byte floats. It then applies the functional equations,
 * doubling the number of points each time.
 */
