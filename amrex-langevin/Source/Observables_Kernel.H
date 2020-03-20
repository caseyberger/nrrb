#ifndef _OBSERVABLES_K_H
#define _OBSERVABLES_K_H

#include "Observables.H"

using namespace amrex;

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real S_tau_Re(int i,int j,int t,int a, Real mu,amrex::Array4<const amrex::Real> const& Lattice){
	Real S_Re = 0.;
	S_Re +=  0.5 * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t,Field(a,C::Re));
	S_Re += -0.5 * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t,Field(a,C::Im));
	S_Re += -0.5 * std::exp(mu) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(a,C::Re));
	S_Re +=  0.5 * std::exp(mu) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(a,C::Im));
	for (int b=1; b<=2; b++){
		S_Re += 0.5 * std::exp(mu) * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(b,C::Im));
		S_Re += 0.5 * std::exp(mu) * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(b,C::Re));
	}
	return S_Re;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real S_tau_Im(int i,int j,int t,int a,Real mu, amrex::Array4<const amrex::Real> const& Lattice){
	Real S_Im = 0.;
	S_Im +=                  Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t,Field(a,C::Im));
	S_Im += -0.5 * std::exp(mu) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(a,C::Im));
	S_Im += -0.5 * std::exp(mu) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(a,C::Re));
	for (int b=1; b<=2; b++){
		S_Im += -0.5 * std::exp(mu) * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(b,C::Re));
		S_Im +=  0.5 * std::exp(mu) * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(b,C::Im));
	}
	return S_Im;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real S_del_Re(int i,int j,int t,int a,Real m, amrex::Array4<const amrex::Real> const& Lattice,const amrex::GeometryData& geom){
	const auto domain_box = geom.Domain();
    const auto domain_lo = amrex::lbound(domain_box);
    const auto domain_hi = amrex::ubound(domain_box);
    Real S_Re = 0.;
	S_Re +=  0.5 * m * (AMREX_SPACEDIM-1) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t,Field(a,C::Re));
	S_Re += -0.5 * m * (AMREX_SPACEDIM-1) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t,Field(a,C::Im));
	for (int d = 1; d <= AMREX_SPACEDIM-1; d++) { //loop over adjacent spatial sites!
		//if r is on the border, this term is zero because it's a derivative of \phi_{a,r}\phi_{a,r+i}^{R}
		// -(1/2m)\phi_{a,i+1}^{R} - (1/2m)\phi_{a,ii1}^{R}
        if (d == 1)
        {
            if (!(i == domain_lo.x || i == domain_hi.x))
            {
				S_Re += -0.25 * m * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i+1,j,t,Field(a,C::Re));
				S_Re +=  0.25 * m * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i+1,j,t,Field(a,C::Im));
				S_Re += -0.25 * m * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i-1,j,t,Field(a,C::Re));
				S_Re +=  0.25 * m * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i-1,j,t,Field(a,C::Im));
				for (int b=1; b<=2; b++){
					S_Re += 0.25 * m * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i+1,j,t,Field(b,C::Im));
					S_Re += 0.25 * m * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i+1,j,t,Field(b,C::Re));
					S_Re += 0.25 * m * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i-1,j,t,Field(b,C::Im));
					S_Re += 0.25 * m * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i-1,j,t,Field(b,C::Re));
				}//loop over b
            }//checking x boundaries
        }//loop over x dim
        else if (d == 2)
        {
            if (!(j == domain_lo.y || j == domain_hi.y))
            {
				S_Re += -0.25 * m * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j+1,t,Field(a,C::Re));
				S_Re +=  0.25 * m * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j+1,t,Field(a,C::Im));
				S_Re += -0.25 * m * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j-1,t,Field(a,C::Re));
				S_Re +=  0.25 * m * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j-1,t,Field(a,C::Im));
				for (int b=1; b<=2; b++){
					S_Re += 0.25 * m * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j+1,t,Field(b,C::Im));
					S_Re += 0.25 * m * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j+1,t,Field(b,C::Re));
					S_Re += 0.25 * m * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j-1,t,Field(b,C::Im));
					S_Re += 0.25 * m * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j-1,t,Field(b,C::Re));
				}//loop over b
            }//checking y boundaries
        }//loop over y dim
	}//derivative loop
	return S_Re;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real S_del_Im(int i,int j,int t,int a, Real m, amrex::Array4<const amrex::Real> const& Lattice,const amrex::GeometryData& geom){
	const auto domain_box = geom.Domain();
    const auto domain_lo = amrex::lbound(domain_box);
    const auto domain_hi = amrex::ubound(domain_box);
    Real S_Im = 0.;
	S_Im += m * (AMREX_SPACEDIM-1) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t,Field(a,C::Im));
	for (int d = 1; d <= AMREX_SPACEDIM-1; d++) { //loop over adjacent spatial sites!
		//if r is on the border, this term is zero because it's a derivative of \phi_{a,r}\phi_{a,r+i}^{R}
		// -(1/2m)\phi_{a,i+1}^{R} - (1/2m)\phi_{a,ii1}^{R}
		if (d == 1)
        {
            if (!(i == domain_lo.x || i == domain_hi.x))
            {
 				S_Im += -0.25 * m * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i+1,j,t,Field(a,C::Im));
				S_Im += -0.25 * m * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i+1,j,t,Field(a,C::Re));
				S_Im += -0.25 * m * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i-1,j,t,Field(a,C::Im));
				S_Im += -0.25 * m * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i-1,j,t,Field(a,C::Re));
				for (int b=1; b<=2; b++){
					S_Im += -0.25 * m * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i+1,j,t,Field(b,C::Re));
					S_Im +=  0.25 * m * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i+1,j,t,Field(b,C::Im));
					S_Im += -0.25 * m * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i-1,j,t,Field(b,C::Re));
					S_Im +=  0.25 * m * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i-1,j,t,Field(b,C::Im));
				}//loop over b
            }//checking x boundaries
        }//loop over x dim
        else if (d == 2)
        {
            if (!(j == domain_lo.y || j == domain_hi.y))
            {
				S_Im += -0.25 * m * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j+1,t,Field(a,C::Im));
				S_Im += -0.25 * m * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j+1,t,Field(a,C::Re));
				S_Im += -0.25 * m * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j-1,t,Field(a,C::Im));
				S_Im += -0.25 * m * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j-1,t,Field(a,C::Re));
				for (int b=1; b<=2; b++){
					S_Im += -0.25 * m * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j+1,t,Field(b,C::Re));
					S_Im +=  0.25 * m * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j+1,t,Field(b,C::Im));
					S_Im += -0.25 * m * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j-1,t,Field(b,C::Re));
					S_Im +=  0.25 * m * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j-1,t,Field(b,C::Im));
				}//loop over b
            }//checking y boundaries
        }//loop over y dim
	}//derivative loop
	return S_Im;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real S_trap_Re(int i,int j,int t,int a, Real w_t,const Real r2, amrex::Array4<const amrex::Real> const& Lattice){
	Real S_Re = 0.;
	S_Re +=  0.25 * w_t * w_t * r2 * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(a,C::Re));
	S_Re += -0.25 * w_t * w_t * r2 * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(a,C::Im));
	for (int b=1; b<=2; b++){
		S_Re += -0.25 * w_t * w_t * r2 * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(b,C::Im));
		S_Re += -0.25 * w_t * w_t * r2 * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(b,C::Re));
	}//loop over b
	return S_Re;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real S_trap_Im(int i,int j,int t,int a,Real w_t,const Real r2, amrex::Array4<const amrex::Real> const& Lattice){
	Real S_Im = 0.;
	S_Im += 0.25 * w_t * w_t * r2 * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(a,C::Im));
	S_Im += 0.25 * w_t * w_t * r2 * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(a,C::Re));
	for (int b=1; b<=2; b++){
		S_Im +=  0.25 * w_t * w_t * r2 * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(b,C::Re));
		S_Im += -0.25 * w_t * w_t * r2 * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(b,C::Im));
	}//loop over b
	return S_Im;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real S_w_Re(int i,int j,int t,int a,Real w,const Real x,const Real y, amrex::Array4<const amrex::Real> const& Lattice){
	Real S_Re = 0.;
	for (int b = 1; b<=2; b++){
		S_Re +=  0.5 * w * epsilon(a,b) * (x - y) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(b,C::Re));
		S_Re += -0.5 * w * epsilon(a,b) * (x - y) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(b,C::Im));
		S_Re += -0.5 * w * epsilon(a,b) * x * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j-1,t-1,Field(b,C::Re));
		S_Re +=  0.5 * w * epsilon(a,b) * x * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j-1,t-1,Field(b,C::Im));
		S_Re +=  0.5 * w * epsilon(a,b) * y * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i-1,j,t-1,Field(b,C::Re));
		S_Re += -0.5 * w * epsilon(a,b) * y * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i-1,j,t-1,Field(b,C::Im));
		S_Re +=  0.5 * w * (x - y) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(b,C::Im));
		S_Re +=  0.5 * w * (x - y) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(b,C::Re));
		S_Re += -0.5 * w * x * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j-1,t-1,Field(b,C::Im));
		S_Re += -0.5 * w * x * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j-1,t-1,Field(b,C::Re));
		S_Re +=  0.5 * w * y * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i-1,j,t-1,Field(b,C::Im));
		S_Re +=  0.5 * w * y * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i-1,j,t-1,Field(b,C::Re));
	}

	return S_Re;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real S_w_Im(int i,int j,int t,int a,Real w,const Real x,const Real y, amrex::Array4<const amrex::Real> const& Lattice){
	Real S_Im = 0.;
	for (int b = 1; b<=2; b++){
		S_Im +=  0.5 * w * epsilon(a,b) * (x-y) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(b,C::Im));
		S_Im +=  0.5 * w * epsilon(a,b) * (x-y) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(b,C::Re));
		S_Im += -0.5 * w * epsilon(a,b) * x * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j-1,t-1,Field(b,C::Im));
		S_Im += -0.5 * w * epsilon(a,b) * x * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j-1,t-1,Field(b,C::Re));
		S_Im +=  0.5 * w * epsilon(a,b) * y * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i-1,j,t-1,Field(b,C::Im));
		S_Im +=  0.5 * w * epsilon(a,b) * y * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i-1,j,t-1,Field(b,C::Re));
		S_Im +=  0.5 * w * (y-x) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(b,C::Re));
		S_Im += -0.5 * w * (y-x) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(b,C::Im));
		S_Im +=  0.5 * w * x * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j-1,t-1,Field(b,C::Re));
		S_Im += -0.5 * w * x * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j-1,t-1,Field(b,C::Im));
		S_Im += -0.5 * w * y * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i-1,j,t-1,Field(b,C::Re));
		S_Im +=  0.5 * w * y * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i-1,j,t-1,Field(b,C::Im));
	}
	return S_Im;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real S_int_Re(int i,int j,int t,int a,Real l, amrex::Array4<const amrex::Real> const& Lattice){
	Real S_Re = 0.;
	for (int b=1; b<=2; b++){
		S_Re +=  0.25 * l * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(b,C::Re)) * Lattice(i,j,t-1,Field(b,C::Re));
		S_Re += -0.25 * l * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(b,C::Re)) * Lattice(i,j,t-1,Field(b,C::Re));
		S_Re +=  0.25 * l * ( Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(b,C::Im)) * Lattice(i,j,t-1,Field(b,C::Im)) - Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(b,C::Im)) * Lattice(i,j,t-1,Field(b,C::Im)) );
		S_Re +=  0.5  * l * ( Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(a,C::Im)) * Lattice(i,j,t,Field(b,C::Im)) * Lattice(i,j,t-1,Field(b,C::Im)) - Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(a,C::Im)) * Lattice(i,j,t,Field(b,C::Re)) * Lattice(i,j,t-1,Field(b,C::Re)) );
		S_Re +=  0.5  * l * ( Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(a,C::Re)) * Lattice(i,j,t,Field(b,C::Re)) * Lattice(i,j,t-1,Field(b,C::Re)) - Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(a,C::Re)) * Lattice(i,j,t,Field(b,C::Im)) * Lattice(i,j,t-1,Field(b,C::Im)) );
		S_Re += -0.5  * l * ( Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(a,C::Im)) * Lattice(i,j,t,Field(b,C::Re)) * Lattice(i,j,t-1,Field(b,C::Im)) + Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(a,C::Re)) * Lattice(i,j,t,Field(b,C::Im)) * Lattice(i,j,t-1,Field(b,C::Re)) );
		S_Re += -0.5  * l * ( Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(a,C::Im)) * Lattice(i,j,t,Field(b,C::Im)) * Lattice(i,j,t-1,Field(b,C::Re)) + Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(a,C::Re)) * Lattice(i,j,t,Field(b,C::Re)) * Lattice(i,j,t-1,Field(b,C::Im)) );
		S_Re += l * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(b,C::Re)) * Lattice(i,j,t-1,Field(b,C::Im));
		S_Re += 0.25 * l * epsilon(a,b) * ( Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(a,C::Re)) * Lattice(i,j,t-1,Field(a,C::Re)) * Lattice(i,j,t,Field(b,C::Im)) - Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(a,C::Im)) * Lattice(i,j,t-1,Field(a,C::Im)) * Lattice(i,j,t,Field(b,C::Im)) );
		S_Re += 0.25 * l * epsilon(a,b) * ( Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(a,C::Re)) * Lattice(i,j,t-1,Field(a,C::Re)) * Lattice(i,j,t,Field(b,C::Re)) - Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(a,C::Im)) * Lattice(i,j,t-1,Field(a,C::Im)) * Lattice(i,j,t,Field(b,C::Re)) );
		S_Re += 0.25 * l * epsilon(a,b) * ( Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(a,C::Re)) * Lattice(i,j,t-1,Field(b,C::Im)) - Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(a,C::Re)) * Lattice(i,j,t-1,Field(b,C::Im)) );
		S_Re += 0.25 * l * epsilon(a,b) * ( Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(a,C::Im)) * Lattice(i,j,t-1,Field(b,C::Re)) - Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(a,C::Im)) * Lattice(i,j,t-1,Field(b,C::Re)) );
		S_Re += 0.5  * l * epsilon(a,b) * ( Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(a,C::Im)) * Lattice(i,j,t-1,Field(a,C::Re)) * Lattice(i,j,t,Field(b,C::Re)) - Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(a,C::Im)) * Lattice(i,j,t-1,Field(a,C::Re)) * Lattice(i,j,t,Field(b,C::Im)) );
		S_Re += 0.5  * l * epsilon(a,b) * ( Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(a,C::Re)) * Lattice(i,j,t,Field(b,C::Im)) * Lattice(i,j,t-1,Field(b,C::Im)) - Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(a,C::Re)) * Lattice(i,j,t,Field(b,C::Re)) * Lattice(i,j,t-1,Field(b,C::Re)) );
	}//loop over b
	return S_Re;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real S_int_Im(int i,int j,int t,int a,Real l, amrex::Array4<const amrex::Real> const& Lattice){
	Real S_Im = 0.;
	for (int b=1; b<=2; b++){
		//S_int
		S_Im += 0.5 * l * ( Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(b,C::Im)) * Lattice(i,j,t-1,Field(b,C::Im)) - Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(b,C::Re)) * Lattice(i,j,t-1,Field(b,C::Re)) );
		S_Im += 0.5 * l * ( Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(b,C::Re)) * Lattice(i,j,t-1,Field(b,C::Im)) - Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(b,C::Re)) * Lattice(i,j,t-1,Field(b,C::Im)) );
		S_Im += 0.5 * l * ( Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(a,C::Re)) * Lattice(i,j,t,Field(b,C::Re)) * Lattice(i,j,t-1,Field(b,C::Re)) - Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(a,C::Im)) * Lattice(i,j,t,Field(b,C::Im)) * Lattice(i,j,t-1,Field(b,C::Im)) );
		S_Im += 0.5 * l * ( Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(a,C::Re)) * Lattice(i,j,t,Field(b,C::Re)) * Lattice(i,j,t-1,Field(b,C::Im)) - Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(a,C::Im)) * Lattice(i,j,t,Field(b,C::Im)) * Lattice(i,j,t-1,Field(b,C::Re)) );
		S_Im += 0.5 * l * ( Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(a,C::Im)) * Lattice(i,j,t,Field(b,C::Re)) * Lattice(i,j,t-1,Field(b,C::Re)) - Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(a,C::Re)) * Lattice(i,j,t,Field(b,C::Im)) * Lattice(i,j,t-1,Field(b,C::Im)) );
		S_Im += 0.5 * l * ( Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(a,C::Re)) * Lattice(i,j,t,Field(b,C::Im)) * Lattice(i,j,t-1,Field(b,C::Re)) - Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(a,C::Im)) * Lattice(i,j,t,Field(b,C::Re)) * Lattice(i,j,t-1,Field(b,C::Im)) );
		S_Im += 0.25 * l * epsilon(a,b) * ( Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(a,C::Re)) * Lattice(i,j,t-1,Field(b,C::Re)) - Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(a,C::Im)) * Lattice(i,j,t-1,Field(b,C::Im)) );
		S_Im += 0.25 * l * epsilon(a,b) * ( Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(a,C::Im)) * Lattice(i,j,t-1,Field(a,C::Im)) * Lattice(i,j,t,Field(b,C::Re)) - Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(a,C::Re)) * Lattice(i,j,t-1,Field(a,C::Re)) * Lattice(i,j,t,Field(b,C::Re)) );
		S_Im += 0.25 * l * epsilon(a,b) * ( Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(a,C::Re)) * Lattice(i,j,t-1,Field(a,C::Re)) * Lattice(i,j,t,Field(b,C::Im)) - Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(a,C::Im)) * Lattice(i,j,t-1,Field(a,C::Im)) * Lattice(i,j,t,Field(b,C::Im)) );
		S_Im += 0.25 * l * epsilon(a,b) * ( Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(a,C::Im)) * Lattice(i,j,t-1,Field(b,C::Im)) - Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(a,C::Re)) * Lattice(i,j,t-1,Field(b,C::Re)) );
		S_Im +=  0.5 * l * epsilon(a,b) * ( Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(a,C::Im)) * Lattice(i,j,t-1,Field(a,C::Re)) * Lattice(i,j,t,Field(b,C::Re)) + Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(a,C::Im)) * Lattice(i,j,t-1,Field(a,C::Re)) * Lattice(i,j,t,Field(b,C::Im)) );
		S_Im += -0.5 * l * epsilon(a,b) * ( Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(a,C::Re)) * Lattice(i,j,t-1,Field(b,C::Im)) + Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(a,C::Im)) * Lattice(i,j,t-1,Field(b,C::Re)) );
	}//loop over b
	return S_Im;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real Theta(amrex::Array4<const amrex::Real> const& Lattice,
           int i, int j, int k)
{
    // tan(theta) = (phi_1_Im + phi_2_Re)/(phi_1_Re - phi_2_Im);
    Real theta = std::atan2(Lattice(i,j,k,Field(1,C::Im)) + Lattice(i,j,k,Field(2,C::Re)),
                            Lattice(i,j,k,Field(1,C::Re)) - Lattice(i,j,k,Field(2,C::Im)));

    // return theta in the range [0, 2*pi]
    if (theta < 0.0)
        theta = theta + Constants::TwoPi;

    return theta;
};

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real DeltaTheta(amrex::Array4<const amrex::Real> const& Lattice,
                int ilp, int jlp, int klp,
                int il, int jl, int kl)
{
    // Returns Theta(ilp, jlp, klp) - Theta(il, jl, kl) = theta_l+1 - theta_l
    Real dtheta = Theta(Lattice, ilp, jlp, klp) - Theta(Lattice, il, jl, kl);

    // Put dtheta into the range [-pi, 0] or [0, pi]
    if (dtheta < -Constants::Pi)
        dtheta +=  Constants::TwoPi;
    else if (dtheta > Constants::Pi)
        dtheta += -Constants::TwoPi;

    return dtheta;
};

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real Circulation(const int i, const int j, const int t,
                 amrex::Array4<const amrex::Real> const& Lattice,
				 const amrex::GeometryData& geom, int radius) {
	// Domain geometry
    const auto domain_box = geom.Domain();
    const auto domain_lo = amrex::lbound(domain_box);
    const auto domain_hi = amrex::ubound(domain_box);

	// For this to work properly, the following must be true:
	// - the domain sizes in x and y must be odd
	// - the center of the domain must satisfy the following:
    const int x_center = 0.5 * (domain_lo.x + domain_hi.x);
    const int y_center = 0.5 * (domain_lo.y + domain_hi.y);

	// Loop left/right edges in the x-dimension
	const int i_left = x_center - radius;
    const int i_right = x_center + radius;

	// Loop bottom/top edges in the y-dimension
    const int j_bottom = y_center - radius;
    const int j_top = y_center + radius;

	// Number of points in the t-dimension in the domain
	const long Nt = domain_box.length(AMREX_SPACEDIM-1);

	// Return 0.0 for this box if t is not in the center of the domain
	if (t != Nt/2)
	{
		return 0.0;
	}

	// Initialize circulation sum for this cell
	Real circ_sum = 0.0;

	// We are summing contributions from theta_t_l+1 - theta_t_l
	// where l denotes (x, y) for a loop lattice site
	// and l+1 denotes (x, y) for the next lattice site on the loop
	// chosen by a loop traversal in the -z direction using the right-hand-rule.
	//
	// The loop starts from the point S, around the center of the domain C:
	//
	// ^ y
	// |      ^------->
	// |      |       |
	// |      |   C   |
	// |      |       |
	// |      S<------v
	// |
	// .-----------> x

	// This left/right section includes the 4 corners of the loop
    // if at left and left cell at this y is within the box, add its phase difference
    if (i == i_left && i_left >= lo.x && i_left <= hi.x) {
        if (j >= j_bottom && j < j_top) {
            circ_sum += DeltaTheta(Lattice, i_left, j+1, t, i_left, j, t);
        }
        else if (j == j_top) {
            circ_sum += DeltaTheta(Lattice, i_left+1, j, t, i_left, j, t);
        }
    }

    // if at right and right cell at this y is within the box, add its phase difference
    if (i == i_right && i_right >= lo.x && i_right <= hi.x) {
        if (j > j_bottom && j <= j_top) {
            circ_sum += DeltaTheta(Lattice, i_right, j-1, t, i_right, j, t);
        }
        else if (j == j_bottom) {
            circ_sum += DeltaTheta(Lattice, i_right-1, j, t, i_right, j, t);
        }
    }

	// Add contributions from bottom/top edges of the loop
    // In the following `if` statement, we do not include corners
    // (i.e. we use `>` and `<` instead of `>=` and `<=`)
    // because the above left/right sections already accounted for them.
    if (i > i_left && i < i_right) {
        // if top cell at this x is within the box, add its phase difference
        if (j == j_top && j_top >= lo.y && j_top <= hi.y) {
            circ_sum += DeltaTheta(Lattice, i+1, j_top, t, i, j_top, t);
        }

        // if bottom cell at this x is within the box, add its phase difference
        if (j == j_bottom && j_bottom >= lo.y && j_bottom <= hi.y) {
            circ_sum += DeltaTheta(Lattice, i-1, j_bottom, t, i, j_bottom, t);
        }
    }

	// Return the circulation contribution from this cell
	Real circ = circ_sum / Constants::TwoPi;
	return circ;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::GpuArray<amrex::Real, Obs::NumObservables>
local_observables(const int i, const int j, const int t,
                  const amrex::Array4<const amrex::Real>& Lattice,
                  const amrex::GeometryData& geom,
                  const NRRBParameters nrrb_parm) {
	Real dens_Re = 0;
	Real dens_Im = 0;
	Real phisq_Re = 0.;
	Real phisq_Im = 0.;
	Real Lz_Re = 0.;
	Real Lz_Im = 0.;
    Real S_Re = 0.;
    Real S_Im = 0.;

	//modify parameters by dtau
	const Real mu = nrrb_parm.dtau * nrrb_parm.mu;
	const Real m = nrrb_parm.m / nrrb_parm.dtau;
	const Real w_t = nrrb_parm.dtau * nrrb_parm.w_t;
	const Real w = nrrb_parm.dtau * nrrb_parm.w;
	const Real l = nrrb_parm.dtau * nrrb_parm.l;

	//define spacing parameters (x, y, r^2, etc)
	const auto domain_xlo = geom.ProbLo();
    const auto domain_xhi = geom.ProbHi();

    const auto domain_box = geom.Domain();
    const auto domain_lo = amrex::lbound(domain_box);
    const auto domain_hi = amrex::ubound(domain_box);
	long int domain_volume = domain_box.volume();

    const Real x_center = 0.5 * (domain_xlo[0] + domain_xhi[0]);
    const Real y_center = 0.5 * (domain_xlo[1] + domain_xhi[1]);
    const Real dx_cell = geom.CellSize(0);
    const Real dy_cell = geom.CellSize(1);

    const Real x_cell = domain_xlo[0] + dx_cell * (i + 0.5 - domain_lo.x);
    const Real y_cell = domain_xlo[1] + dy_cell * (j + 0.5 - domain_lo.y);
    const Real r2 = std::pow(x_cell - x_center, 2) + std::pow(y_cell - y_center, 2);

    for (int a = 1; a <= 2; a++){
        //Field modulus squared
        phisq_Re +=  0.5 * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t,Field(a,C::Re));
        phisq_Re += -0.5 * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t,Field(a,C::Im));
        phisq_Im += Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t,Field(a,C::Im));
        for (int b = 1; b <= 2; b++){
            //Density
            dens_Re +=  delta(a,b)   * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(b,C::Re));
            dens_Re += -delta(a,b)   * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(b,C::Im));
            dens_Re += -epsilon(a,b) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(b,C::Re));
            dens_Re += -epsilon(a,b) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(b,C::Im));
            dens_Im +=  delta(a,b)   * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(b,C::Re));
            dens_Im +=  delta(a,b)   * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(b,C::Im));
            dens_Im +=  epsilon(a,b) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(b,C::Re));
            dens_Im += -epsilon(a,b) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(b,C::Im));
        }//loop over b = 1,2
    }//loop over a = 1,2

    //Angular momentum
    if (AMREX_SPACEDIM == 3){
        //define x and y
        const Real x = x_cell - x_center;
        const Real y = y_cell - y_center;
        for (int a = 1; a <= 2; a++){//loop over a
            //real sum over a only
            Lz_Re +=  0.5 * x * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j-1,t-1,Field(a,C::Im));
            Lz_Re +=  0.5 * x * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j-1,t-1,Field(a,C::Re));
            Lz_Re += -0.5 * y * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i-1,j,t-1,Field(a,C::Im));
            Lz_Re += -0.5 * y * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i-1,j,t-1,Field(a,C::Re));
            Lz_Re +=        y * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t,Field(a,C::Im));
            Lz_Re +=       -x * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t,Field(a,C::Im));
            //imaginary sum over a only
            Lz_Im += -0.5 * x * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j-1,t-1,Field(a,C::Re));
            Lz_Im +=  0.5 * x * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j-1,t-1,Field(a,C::Im));
            Lz_Im +=  0.5 * y * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i-1,j,t-1,Field(a,C::Re));
            Lz_Im += -0.5 * y * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i-1,j,t-1,Field(a,C::Im));
            Lz_Im += -0.5 * y * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t,Field(a,C::Re));
            Lz_Im +=  0.5 * y * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t,Field(a,C::Im));
            Lz_Im +=  0.5 * x * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t,Field(a,C::Re));
            Lz_Im += -0.5 * x * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t,Field(a,C::Im));
            for (int b = 1; b <= 2; b++){//loop over b
                //real sum over a and b
                Lz_Re +=  0.5 * epsilon(a,b) * x * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j-1,t-1,Field(b,C::Re));
                Lz_Re += -0.5 * epsilon(a,b) * x * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j-1,t-1,Field(b,C::Im));
                Lz_Re +=  0.5 * epsilon(a,b) * y * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i-1,j,t-1,Field(b,C::Im));
                Lz_Re += -0.5 * epsilon(a,b) * y * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i-1,j,t-1,Field(b,C::Re));
                //imaginary sum over a and b
                Lz_Im += -0.5 * epsilon(a,b) * x * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j-1,t-1,Field(b,C::Im));
                Lz_Im += -0.5 * epsilon(a,b) * x * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j-1,t-1,Field(b,C::Re));
                Lz_Im +=  0.5 * epsilon(a,b) * y * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i-1,j,t-1,Field(b,C::Im));
                Lz_Im +=  0.5 * epsilon(a,b) * y * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i-1,j,t-1,Field(b,C::Re));
            }//loop over b for Lz
        }//loop over a for Lz
    }//if dim+1 = 3

    //ACTION
    for (int a = 1; a <= 2; a++){
        //S_tau
        S_Re += S_tau_Re(i,j,t,a,mu,Lattice);
        S_Im += S_tau_Im(i,j,t,a,mu,Lattice);

        //S_del
        S_Re += S_del_Re(i,j,t,a,m,Lattice,geom);
        S_Im += S_del_Im(i,j,t,a,m,Lattice,geom);

        //S_trap
        S_Re += S_trap_Re(i,j,t,a,w_t,r2,Lattice);
        S_Im += S_trap_Im(i,j,t,a,w_t,r2,Lattice);

        //S_w
        if (AMREX_SPACEDIM == 3){
            const Real x = x_cell - x_center;
            const Real y = y_cell - y_center;
            S_Re += S_w_Re(i,j,t,a,w,x,y,Lattice);
            S_Im += S_w_Im(i,j,t,a,w,x,y,Lattice);
        }

        //S_int
        S_Re += S_int_Re(i,j,t,a,l,Lattice);
        S_Im += S_int_Im(i,j,t,a,l,Lattice);
    }//loop over a (for the Action)

    dens_Re = 0.5 * std::exp(mu) * dens_Re;
    dens_Im = 0.5 * std::exp(mu) * dens_Im;

    GpuArray<Real, Obs::NumObservables> observables;

    observables[Obs::PhiSqRe] = phisq_Re / domain_volume;
    observables[Obs::PhiSqIm] = phisq_Im / domain_volume;
    observables[Obs::DensRe]  = dens_Re  / domain_volume;
    observables[Obs::DensIm]  = dens_Im  / domain_volume;
    observables[Obs::LzRe]    = Lz_Re    / domain_volume;
    observables[Obs::LzIm]    = Lz_Im    / domain_volume;
    observables[Obs::SRe]     = S_Re     / domain_volume;
    observables[Obs::SIm]     = S_Im     / domain_volume;

    // compute circulation contributions from this cell
    observables[Obs::Circ1] = Circulation(i, j, t, Lattice, geom, nrrb_parm.circulation_radius_1);
    observables[Obs::Circ2] = Circulation(i, j, t, Lattice, geom, nrrb_parm.circulation_radius_2);

    return observables;
}

#endif
