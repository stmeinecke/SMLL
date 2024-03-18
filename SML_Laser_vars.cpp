class vars{ 
	 public:
		std::complex<double> E;	
    double rho;
    double N;

		vars(){this->setTo(0.0);}

		void add(vars *v, vars *r){
			for(int i=0; i<sizeof(vars)/sizeof(double); i+=1){
				((double*)r)[i] = ((double*)this)[i] + ((double*)v)[i];
			}
		} 

		void mult(vars *r, double m){
			for(int i=0; i<sizeof(vars)/sizeof(double); i+=1){
				((double*)r)[i] = ((double*)this)[i]*m ;
			}
		}

		void setTo(double s){
			for(int i=0; i<sizeof(vars)/sizeof(double); i+=1){
				((double*)this)[i] = s;
			}
		}

};

namespace IND{
  
  int E = 0;
  int rho = 2;
  int N = 3;
  
}
