#include<iostream>
#include<utility>
#include<vector>
#include<cmath>
#include<omp.h>

#include<boost/multiprecision/cpp_int.hpp>
//#include<boost/multiprecision/gmp.hpp>

uint_fast16_t const constexpr prime[256]={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541,547,557,563,569,571,577,587,593,599,601,607,613,617,619,631,641,643,647,653,659,661,673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,773,787,797,809,811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,941,947,953,967,971,977,983,991,997,1009,1013,1019,1021,1031,1033,1039,1049,1051,1061,1063,1069,1087,1091,1093,1097,1103,1109,1117,1123,1129,1151,1153,1163,1171,1181,1187,1193,1201,1213,1217,1223,1229,1231,1237,1249,1259,1277,1279,1283,1289,1291,1297,1301,1303,1307,1319,1321,1327,1361,1367,1373,1381,1399,1409,1423,1427,1429,1433,1439,1447,1451,1453,1459,1471,1481,1483,1487,1489,1493,1499,1511,1523,1531,1543,1549,1553,1559,1567,1571,1579,1583,1597,1601,1607,1609,1613,1619};

float const constexpr primelog[256]={1.000,1.585,2.322,2.807,3.459,3.700,4.087,4.248,4.524,4.858,4.954,5.209,5.358,5.426,5.55,5.728,5.883,5.931,6.066,6.150,6.190,6.304,6.375,6.476,6.600,6.658,6.687,6.741,6.768,6.20,6.989,7.033,7.098,7.119,7.219,7.238,7.295,7.349,7.384,7.435,7.484,7.500,7.577,7.592,7622,7.637,7.721,7.801,7.827,7.839,7.864,7.901,7.913,7.972,8.006,8.039,8.071,8.082,8.114,.134,8.145,8.195,8.262,8.281,8.290,8.308,8.371,8.397,8.439,8.447,8.464,8.488,8.520,8.543,8.566,8.581,8.604,8.633,8.647,8.676,8.711,8.718,8.752,8.758,8.778,8.791,8.811,8.836,8.89,8.855,8.867,8.904,8.928,8.940,8.963,8.974,8.992,9.025,9.031,9.079,9.095,9.122,9.137,9.52,9.157,9.172,9.197,9.212,9.226,9.231,9.246,9.260,9.269,9.274,9.301,9.324,9.329,9.338,9351,9.364,9.369,9.394,9.403,9.416,9.433,9.453,9.470,9.490,9.506,9.518,9.529,9.537,9.553,.564,9.572,9.587,9.594,9.620,9.638,9.660,9.664,9.681,9.685,9.692,9.695,9.713,9.736,9.743,9.747,9.753,9.776,9.783,9.786,9.793,9.825,9.831,9.844,9.860,9.872,9.878,9.887,9.896,9.97,9.923,9.932,9.941,9.953,9.961,9.979,9.984,9.993,9.996,10.010,10.013,10.021,10.035,10.08,10.051,10.054,10.062,10.086,10.091,10.094,10.099,10.107,10.115,10.125,10.133,10.141,10169,10.171,10.184,10.194,10.206,10.213,10.220,10.230,10.244,10.249,10.256,10.263,10.266,0.273,10.287,10.298,10.319,10.321,10.325,10.332,10.334,10.341,10.345,10.348,10.352,10.365,10.367,10.374,10.410,10.417,10.423,10.431,10.450,10.460,10.475,10.479,10.481,10.485,10.491,10.499,10.53,10.505,10.511,10.523,10.532,10.534,10.538,10.540,10.544,10.550,10.561,10.573,10.580,10.92,10.597,10.601,10.606,10.614,10.617,10.625,10.628,10.641,10.645,10.650,10.652,10.656,10.661};

int_fast32_t const N=64;

using intA = boost::multiprecision::cpp_int;

uint_fast16_t inline depth(std::vector<std::pair<uint_fast64_t,uint_fast64_t>> pc, uint_fast64_t base_){
	uint_fast16_t d=1;

	intA n=1;
	
	for(uint_fast32_t i=0;i<pc.size();i++){ intA ni=pc[i].first; n=n*boost::multiprecision::pow(ni,pc[i].second); }
	
	intA nit=n;
	intA const base=base_;
	do{
		intA product=1;
	
		while(nit!=0 && product!=0 && product%base!=0){
		       	intA digit=nit%base;
		
			if(digit==0){ product=0; break;}
			else if(digit!=1) product=product*digit; 
		
			nit=nit/base; 
		}
	
		nit=product;
		d++;
	}while(nit/base!=0);
	return d;
}

int main(int argc, char* argv[]){
	uint_fast32_t const Nthreads=omp_get_max_threads();

	omp_set_dynamic(0);
	omp_set_num_threads(Nthreads);

	for(uint_fast64_t base=10;base<11;base++){
		uint_fast16_t max=0;
		std::vector<std::pair<uint_fast64_t,uint_fast64_t>> pc_;
	
		for(uint_fast16_t i=0;i<256;i++){
			if(prime[i]<base) pc_.push_back({prime[i],0});
			else break;
		}

	
		std::vector<std::vector<std::pair<uint_fast64_t,uint_fast64_t>>> pc;
		for(uint_fast32_t i=0;i<Nthreads;i++) pc.push_back(pc_);

		uint_fast32_t total=1;
		for(uint_fast32_t i=0; i<pc[0].size();i++) total=total*N;

		std::vector<std::vector<std::pair<uint_fast64_t,uint_fast64_t>>> list;

		#pragma omp parallel for 
		for(uint_fast32_t i=0;i<total;i++){
			uint_fast32_t tid=omp_get_thread_num();
			uint_fast32_t tmp;

			uint_fast16_t d=depth(pc[tid],base);

			if(d>=max){
				#pragma omp critical
				{
					if(d>max){ list.clear();max=d;}
					if(d>=max) list.push_back(pc[tid]);	
				}
			}

			tmp=i;
			for(uint_fast32_t j=0;j<pc[tid].size();j++){
				pc[tid][j].second=tmp%N;
				tmp=tmp/N;
			}	
		}
		std::cerr << base << '\t' << list.size() << '\t' << max << std::endl;
		for(uint_fast16_t i=0;i<list.size();i++){
			for(uint_fast16_t j=0;j<list[i].size()-1;j++) std::cerr << list[i][j].first << ' ' << list[i][j].second << ", ";
			std::cerr << list[i][list[i].size()-1].first << ' ' << list[i][list[i].size()-1].second << std::endl;
		}
		std::cerr << std::endl;
		
		std::cout << base << '\t' << list.size() << '\t' << max << std::endl;
		for(uint_fast16_t i=0;i<list.size();i++){
			for(uint_fast16_t j=0;j<list[i].size()-1;j++) std::cout << list[i][j].first << ' ' << list[i][j].second << ", ";
			std::cout << list[i][list[i].size()-1].first << ' ' << list[i][list[i].size()-1].second << std::endl;
		}
		std::cout << std::endl;
	}

	return 0;
}
