 
#include "helper.hpp"

/**
 * Determine if i and j are close within a neglible amount i.e.
 * - absolute value of the difference between i and j is smaller than the minimum value of a double
 * - absolute value of the scaled machine epsilon to the units of least precision
 * 
 * The intendend usage in this application is in comparison to 0
 * 
 * @param i         First number to compare
 * @param j         Second number to compare
 * @param ulp       Units of least precision
 * 
 * @return          bool
 */
bool Helper::almost_equal(double i, double j, int ulp) {

    return  fabs(i - j) <= (numeric_limits<double>::epsilon() * fabs(i + j) * ulp)
             || (fabs(i - j) < numeric_limits<double>::min());
   
}

/**
 * Determine if value is inclusively between min and max
 * 
 * @param value     Value to comapre to determine between min and max
 * @param min       Start of valid range for value
 * @param max       End of valid range for value
 * 
 * @return          bool true if value between min and max
 *                  false otherwise
 */
bool Helper::is_between(double value, double min, double max) {
    return (value >= min) && (value <= max);
}

/**
 * Determine if file_name exists on disk
 * 
 * @param file_name     Char array containing dirctory, file and extension to check
 * 
 * @return              bool true if file_name exists on the local disk
 */
bool Helper::is_file_exist(const char *file_name) {
    ifstream infile(file_name);
    return infile.good();
}

/**
 * Safe multiplication of a and b
 * 
 * @param a     first value
 * @param b     second value
 * @return      double a * b
 *              throws invalid_argument in case of overflow
 */
double Helper::times(double a, double b) {
    
    // Note there are precision issues cause the division to be unequal. Can implement a better check
    double x = a * b;
    if (a != 0 && abs(x/a-b) > 1) {
        throw invalid_argument("Mutiplication overflow a:"+to_string(a)+" b:"+to_string(b)+" x:"+to_string(x)+" x/a:"+to_string(x/a));   
    }
    return x;
}

/**
 * Safe divide of a and b
 * 
 * @param a     first value
 * @param b     second value
 * @return      double a / b
 *              throws invalid_argument in case of overflow
 */
double Helper::divide(double a, double b) {
    
    if( b == 0)
        throw invalid_argument("Divide by 0");

    // If the denominator is tiny we can overflow, so I think we need this..

    double x = a / b;
    if (abs(x*b-a) > 1)
        throw invalid_argument("Divide overflow a:"+to_string(a)+" b:"+to_string(b)+" x:"+to_string(x)+" x*b:"+to_string(x*b));   
        
    return x;
}

/**
 * Safe addition of a and b
 * 
 * @param a     first value
 * @param b     second value
 * @return      double a + b
 *              throws invalid_argument in case of overflow
 */
double Helper::add(double a, double b) {

    double x = a + b;
    if (a > 0 && b > 0 && x < 0) {
        throw invalid_argument("Addition overflow :"+to_string(a)+" b:"+to_string(b)+" x:"+to_string(x));    
    }
    if (a < 0 && b < 0 && x > 0) {
        throw invalid_argument("Addition underflow :"+to_string(a)+" b:"+to_string(b)+" x:"+to_string(x));  
    }
    return x;
}

/**
 * Safe exponent raised to a
 * 
 * @param a     value
 * @return      double exp(a)
 *              throws invalid_argument in case of overflow
 */
double Helper::exp(double a) {

    double ret = std::exp(a);

    if (errno) 
        throw invalid_argument("Exponent overflow");   

    return ret;   
}

/**
 * Safe raising of power to base 
 * 
 * @param base      base of operation
 * @param power     power of operation
 * @return          double pow(base, power)
 *                  throws invalid_argument in case of overflow
 */
double Helper::pow(double base, double power) {

    double ret = std::pow(base, power);
    if (errno) 
        throw invalid_argument("Powers overflow");   

    return ret; 
}

/**
 * Algorithm to find the GCD of a and b
 * 
 * @param a       one of the number to compute the GCD
 * @param b       another number to compute the GCD 
 * @return        return the integer value of GCD between a and b
 */
long long Helper::gcd(long long a,  long long b)
{
    if (a == 0) // If a = 0 then GCD(a,b)=b and stop
        return b;
    else if (b == 0)     // If b = 0 then GCD(a,b)=a and stop
        return a;
    if (a < b) 
        return gcd(a, b % a); 
    else
        return gcd(b, a % b);
}


/**
 * Function to convert decimal to Rational
 * 
 * @param coeff     the decimal number for which we will compute the rational for
 * @param numer     parameter by reference: to return the numerator of the rational form
 * @param denom     parameter by reference: to return the denominator of the rational form
 */
void Helper::decimalToRational(double coeff, long long &numer, long long &denom)
{
    double intVal = floor(coeff);           // Fetch integer part of the decimal coefficent
    double fracVal = coeff - intVal;        // Compute the fractional part of the coefficent
    const long long precisionVal = 1e18;    // Precision to the fractional part

    // Calculate GCD of integral equivalent of fractional part and precision value
    long long gcdVal    = gcd(round(fracVal * precisionVal), precisionVal);
 
    // Calculate numerator and denominator of the supplied coefficient
    long long coefff_numerator   = round(fracVal * precisionVal) / gcdVal;
    denom = precisionVal / gcdVal; 
    numer = (long long) (intVal * denom) + coefff_numerator;
}

/**
 * Split string based on specific delimiter
 * 
 * @param str       string to processs
 * @param delim     delimiter to split string
 * @return          vector of split strings
 *                  
 */
vector<string> Helper::split(string str, string delim) {

    size_t pos_start = 0;
    size_t pos_end;
    size_t delim_len = delim.length();

    string token;
    vector<string> res;

    while ((pos_end = str.find(delim, pos_start)) != string::npos) {
        token = str.substr(pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back(token);
    }

    res.push_back(str.substr(pos_start));
    return res;
}

/**
 * Trim left empty characters of a string
 * 
 * @param str       string to processs
 * @return          string with left spaces removed
 *                  
 */
string  Helper::ltrim(string str) {
    return regex_replace(str, regex("^\\s+"), string(""));
}
 
/**
 * Trim right empty characters of a string
 * 
 * @param str       string to processs
 * @return          string with right spaces removed
 *                  
 */
string  Helper::rtrim(string s) {
    return regex_replace(s, regex("\\s+$"), string(""));
}

/**
 * Trim left and right empty characters of a string
 * 
 * @param str       string to processs
 * @return          string with right spaces removed
 *                  
 */ 
string  Helper::trim(string s) {
    return  ltrim( rtrim(s));
}