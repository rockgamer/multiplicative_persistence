#include<iostream>
#include<utility>
#include<vector>
#include<cmath>
#include<omp.h>

#include<boost/multiprecision/cpp_int.hpp>
#include<boost/multiprecision/gmp.hpp>

void f(){
	//apenas um comentario para ocupar espaco
}

uint_fast32_t const prime[256]={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541,547,557,563,569,571,577,587,593,599,601,607,613,617,619,631,641,643,647,653,659,661,673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,773,787,797,809,811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,941,947,953,967,971,977,983,991,997,1009,1013,1019,1021,1031,1033,1039,1049,1051,1061,1063,1069,1087,1091,1093,1097,1103,1109,1117,1123,1129,1151,1153,1163,1171,1181,1187,1193,1201,1213,1217,1223,1229,1231,1237,1249,1259,1277,1279,1283,1289,1291,1297,1301,1303,1307,1319,1321,1327,1361,1367,1373,1381,1399,1409,1423,1427,1429,1433,1439,1447,1451,1453,1459,1471,1481,1483,1487,1489,1493,1499,1511,1523,1531,1543,1549,1553,1559,1567,1571,1579,1583,1597,1601,1607,1609,1613,1619};

int_fast32_t const N=64;

using intA = boost::multiprecision::mpz_int;

uint_fast32_t inline depth(std::vector<std::pair<intA,uint_fast32_t>> pc, uint_fast32_t base_){
	uint_fast32_t d=1;

	intA n=1;
	
	for(uint_fast32_t i=0;i<pc.size();i++) n=n*boost::multiprecision::pow(pc[i].first,pc[i].second);
	
	intA nit=n;
	intA base=base_;
	do{
		intA product=1;
	
		while(nit!=0 && product!=0){
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
	omp_set_num_threads(Nthreads);

	for(uint_fast32_t base=3;base<21;base++){
		uint_fast32_t max=0;
		std::vector<std::pair<uint_fast32_t,uint_fast32_t>> multi_table;	
		std::vector<std::vector<std::pair<intA,uint_fast32_t>>> pc;
		std::vector<std::vector<std::pair<intA,uint_fast32_t>>> list;
		
		bool p=true;
		{
			uint_fast32_t base_=base;	
			std::vector<std::pair<intA,uint_fast32_t>> pc_;	
			std::vector<std::pair<uint_fast32_t,uint_fast32_t>> multi_table_;
			uint_fast32_t idx=0;

			for(uint_fast32_t i=0;i<256;i++){
				if(prime[i]<base){ pc_.push_back({prime[i],0}); multi_table_.push_back({prime[i],0});}
				else break;
			}

			for(uint_fast32_t i=0;i<Nthreads;i++) pc.push_back(pc_);

			do{
				if(base_%multi_table_[idx].first==0){ p=false; multi_table_[idx].second++; base_=base_/multi_table_[idx].first;}
				else idx++;
			}while(base_!=0 && idx<multi_table_.size());
	
			for(uint_fast32_t i=0;i<multi_table_.size();i++){
				if(multi_table_[i].second!=0) multi_table.push_back({i,multi_table_[i].second});
			}
		}

		uint_fast32_t total=std::pow(N,pc[0].size());
	
		#pragma omp parallel for 
		for(uint_fast32_t i=0;i<total;i++){
			uint_fast32_t tid=omp_get_thread_num();
			uint_fast32_t idx;
			bool skip=!p;

			if(skip){
				for(uint_fast32_t j=0;j<multi_table.size();j++){
					if(pc[tid][multi_table[j].first].second<multi_table[j].second){ skip=false; break;}
				}
			}

			if(!skip){
				uint_fast32_t d=depth(pc[tid],base);

				if(d>=max){
					#pragma omp critical
					{
						if(d>max){ list.clear(); max=d;}
						if(d>=max) list.push_back(pc[tid]);	
					}
				}
			}

			idx=i;
			for(uint_fast32_t j=0;j<pc[tid].size();j++){
				pc[tid][j].second=idx%N;
				idx=idx/N;
			}	
		}


		std::clog << base << '\t' << list.size() << '\t' << max << std::endl;
		for(uint_fast32_t i=0;i<list.size();i++){
			for(uint_fast32_t j=0;j<list[i].size()-1;j++) std::clog << list[i][j].first << ' ' << list[i][j].second << ", ";
			std::clog << list[i][list[i].size()-1].first << ' ' << list[i][list[i].size()-1].second << std::endl;
		}
		std::clog << std::endl;
		
		std::cout << base << '\t' << list.size() << '\t' << max << std::endl;
		for(uint_fast32_t i=0;i<list.size();i++){
			for(uint_fast32_t j=0;j<list[i].size()-1;j++) std::cout << list[i][j].first << ' ' << list[i][j].second << ", ";
			std::cout << list[i][list[i].size()-1].first << ' ' << list[i][list[i].size()-1].second << std::endl;
		}
		std::cout << std::endl;
	}

	return 0;
}
