#ifndef GPUFIT_GAMMA_VARIATE1D_CUH_INCLUDED
#define GPUFIT_GAMMA_VARIATE1D_CUH_INCLUDED

/* Description of the calculate_gamma_variate1d function
 * ===================================================
 *
 * Simplified Gamma Variate: M T Madsen 1992 Phys. Med. Biol. 37 1597
 * This function calculates the values of one-dimensional gamma variate model functions
 * and their partial derivatives with respect to the model parameters. 
 *
 * This function makes use of the user information data to pass in the 
 * independent variables (X values) corresponding to the data.  The X values
 * must be of type REAL.
 *
 * Note that if no user information is provided, the (X) coordinate of the 
 * first data value is assumed to be (0.0).  In this case, for a fit size of 
 * M data points, the (X) coordinates of the data are simply the corresponding 
 * array index values of the data array, starting from zero.
 *
 * There are three possibilities regarding the X values:
 *
 *   No X values provided: 
 *
 *       If no user information is provided, the (X) coordinate of the 
 *       first data value is assumed to be (0.0).  In this case, for a 
 *       fit size of M data points, the (X) coordinates of the data are 
 *       simply the corresponding array index values of the data array, 
 *       starting from zero.
 *
 *   X values provided for one fit:
 *
 *       If the user_info array contains the X values for one fit, then 
 *       the same X values will be used for all fits.  In this case, the 
 *       size of the user_info array (in bytes) must equal 
 *       sizeof(REAL) * n_points.
 *
 *   Unique X values provided for all fits:
 *
 *       In this case, the user_info array must contain X values for each
 *       fit in the dataset.  In this case, the size of the user_info array 
 *       (in bytes) must equal sizeof(REAL) * n_points * nfits.
 *
 * Parameters:
 *
 * parameters: An input vector of model parameters.
 *             p[0]: peak(signal at tmax)
 *             p[1]: alpha
 *             p[2]: tmax
 *             p[3]: time of arrival
 *             p[4]: baseline
 * n_fits: The number of fits.
 *
 * n_points: The number of data points per fit.
 *
 * value: An output vector of model function values.
 *
 * derivative: An output vector of model function partial derivatives.
 *
 * point_index: The data point index.
 *
 * fit_index: The fit index.
 *
 * chunk_index: The chunk index. Used for indexing of user_info.
 *
 * user_info: An input vector containing user information.
 *
 * user_info_size: The size of user_info in bytes.
 *
 * Calling the calculate_linear1d function
 * =======================================
 *
 * This __device__ function can be only called from a __global__ function or an other
 * __device__ function.
 *
*/

__device__ void calculate_gamma_variate(
	REAL const * parameters,
	int const n_fits,
	int const n_points,
	REAL * value,
	REAL * derivative,
	int const point_index,
	int const fit_index,
	int const chunk_index,
	char * user_info,
	std::size_t const user_info_size)
{
    // indices

    REAL * user_info_float = (REAL*) user_info;
    REAL x = 0;
    if (!user_info_float)
    {
	x = point_index;
    }
    else if (user_info_size / sizeof(REAL) == n_points)
    {
	x = user_info_float[point_index];
    }
    else if (user_info_size / sizeof(REAL) > n_points)
    {
	int const chunk_begin = chunk_index * n_fits * n_points;
	int const fit_begin = fit_index * n_points;
	x = user_info_float[chunk_begin + fit_begin + point_index];
    }

    /* Parameters
     * parameters: An input vector of model parameters.
     *             p[0]: peak(signal at tmax)
     *             p[1]: alpha
     *             p[2]: tmax
     *             p[3]: time of arrival
     *             p[4]: baseline
     *   
     */
    REAL const * p = parameters;

    //value

    if( x - p[3] <= 0 )
	x = p[3];

    //Constrain denominator in tprime to be positive.
    float p2 = p[2] + p[3];

    //Calculate t prime;
    float s0= (p2 - p[3]);
    float tprime = (x-p[3])/s0;

    // log(x): x must be > 0; 
    value[point_index] = p[0] * pow(tprime, p[1]) * exp(p[1] * (1-tprime)) + p[4];

    // derivatives

    REAL * current_derivatives = derivative + point_index;
    // wrt p[0]
    current_derivatives[0 * n_points] =  pow(tprime, p[1]) * exp((p[1] *(p2 - x))/s0);

    // wrt p[1]
    current_derivatives[1 * n_points] = p[0] * exp(p[1] * (1 - tprime)) * pow(tprime, p[1]) * (1 + log(tprime) - tprime);

    //wrt p2
    current_derivatives[2 * n_points] = p[0] * p[1] * (p[3] - x) * exp(p[1] * (1 - tprime)) * (pow(tprime, p[1]-1)-pow(tprime, p[1]))/pow(s0, 2);

    //wrt p[3]
    current_derivatives[3 * n_points] = p[0] * p[1] * (p2 - x) * exp(p[1] * (1 - tprime))  * (pow(tprime, p[1])-pow(tprime, p[1]-1))/pow(s0, 2);

    // wrt p[4]
    current_derivatives[4 * n_points] = 1;
}

#endif
