#ifndef COLLISION
#define COLLISION

// Necessary includes
#include "macros.cu"
#include "collision.cuh"
#include "solver.cu"

// These files are only included to remove squiggly red lines in VS2010
#include "data_types.cuh"
#include "cuda_runtime.h"

// THIS DECLARATION *SHOULD* BE IN solver.cu FOR CLARITY, HOWEVER COMPILER COMPLAINS IF THIS IS THE CASE
__device__ __constant__ DomainConstant domain_constants;

#define POW4(x) x*x*x*x
#define INVERSEPOW(x) {1./x}

__device__ collision collision_functions[9] = { bgk_collision, bgk_guo_collision, bgk_ntpor_collision, bgk_ntpor_guo_collision,
												mrt_collision, mrt_guo_collision, mrt_ntpor_collision, mrt_ntpor_guo_collision,
												bounceback};

__device__ inline double u_square(Node *current_node)
{
	double value = 0;

	#pragma unroll
	for(int d = 0; d<DIM; d++)
	{
		value += (current_node->u[d]*current_node->u[d]);
	}

	return value*1.5;
}

__device__ inline double e_mul_u(Node *current_node, int *i)
{
	double value = 0;

	#pragma unroll
	for(int d = 0; d<DIM; d++)
	{
		value += (domain_constants.e[d][*i]*current_node->u[d]);
	}

	return value*3.;
}

__device__ __noinline__ void bgk_collision(Node *current_node, double *tau)
{
	double f_eq[Q], u_sq, eu;

	u_sq = u_square(current_node);

	for(int i=0;i<Q;i++)
	{
		eu = e_mul_u(current_node, &i);
		f_eq[i] = current_node->rho*domain_constants.omega[i]*(1.0+eu+(0.5*eu*eu)-u_sq);
	}

	if (current_node->c_smag>0) turbulent_viscosity(current_node, f_eq, tau);

	for(int i = 0; i<Q; i++)
	{
		current_node->f[i] = current_node->f[i] - (1.0/(*tau)) * (current_node->f[i]-f_eq[i]);
	}
}

__device__ __noinline__ void bgk_guo_collision(Node *current_node, double *tau)
{
	double f_eq[Q], u_sq, eu, F_coeff[DIM], force_term[Q];
	int d;
	
	#pragma unroll
	for(d = 0; d<DIM; d++)
	{
		current_node->u[d] = current_node->u[d] + (1/(2*current_node->rho))*current_node->F[d];
	}

	u_sq = u_square(current_node);

	for(int i=0;i<Q;i++)
	{
		#pragma unroll
		for(d = 0; d<DIM; d++)
		{
		#if DIM > 2
			F_coeff[d] = domain_constants.omega[i]*(1-(1/(2*(*tau))))*(((domain_constants.e[d][i]-current_node->u[d])*3)+(domain_constants.e[d][i]*9*((domain_constants.e[0][i]*current_node->u[0])+(domain_constants.e[1][i]*current_node->u[1])+(domain_constants.e[2][i]*current_node->u[2]))));
		#else
			F_coeff[d] = domain_constants.omega[i]*(1-(1/(2*(*tau))))*(((domain_constants.e[d][i]-current_node->u[d])*3)+(domain_constants.e[d][i]*9*((domain_constants.e[0][i]*current_node->u[0])+(domain_constants.e[1][i]*current_node->u[1]))));
		#endif
		}

		force_term[i] = 0;
		#pragma unroll
		for(d = 0; d<DIM; d++)
		{
			force_term[i] += F_coeff[d]*current_node->F[d];
		}

		eu = e_mul_u(current_node, &i);
		f_eq[i] = current_node->rho*domain_constants.omega[i]*(1.0+eu+(0.5*eu*eu)-u_sq);
	}

	if (current_node->c_smag>0) turbulent_viscosity(current_node, f_eq, tau);

	for(int i=0;i<Q;i++)
	{
		current_node->f[i] = current_node->f[i] - (1.0/(*tau)) * (current_node->f[i]-f_eq[i]) + force_term[i];
	}
}

__device__ __noinline__ void bgk_ntpor_collision(Node *current_node, double *tau)
{
	double f_eq[Q], u_sq, eu, collision_bgk, collision_s, tmp[Q];

	u_sq = u_square(current_node);
	for(int i=0;i<Q;i++)
	{
		eu = e_mul_u(current_node, &i);
		f_eq[i] = current_node->rho*domain_constants.omega[i]*(1.0+eu+(0.5*eu*eu)-u_sq);
	}

	if (current_node->c_smag>0) turbulent_viscosity(current_node, f_eq, tau);
	
	for(int i =0;i<Q;i++)
	{
		collision_bgk = (1.0/(*tau)) * (current_node->f[i]-f_eq[i]);
		collision_s = current_node->f[domain_constants.opp[i]]-current_node->f[i];
		tmp[i] = current_node->f[i] - (1-(current_node->B))*collision_bgk + (current_node->B)*collision_s;
	}

	for(int i =0;i<Q;i++)
	{
		current_node->f[i] = tmp[i];
	}

}

__device__ void bgk_ntpor_guo_collision(Node *current_node, double *tau)
{
	double f_eq[Q], u_sq, eu, collision_bgk, collision_s, F_coeff[DIM], force_term[Q], tmp[Q];
	int d;

	#pragma unroll
	for(d = 0; d<DIM; d++)
	{
		current_node->u[d] = current_node->u[d] + (1/(2*current_node->rho))*current_node->F[d];
	}

	u_sq = u_square(current_node);

	for(int i=0;i<Q;i++)
	{
		#pragma unroll
		for(d = 0; d<DIM; d++)
		{
		#if DIM > 2
			F_coeff[d] = domain_constants.omega[i]*(1-(1/(2*(*tau))))*(((domain_constants.e[d][i]-current_node->u[d])*3)+(domain_constants.e[d][i]*9*((domain_constants.e[0][i]*current_node->u[0])+(domain_constants.e[1][i]*current_node->u[1])+(domain_constants.e[2][i]*current_node->u[2]))));
		#else
			F_coeff[d] = domain_constants.omega[i]*(1-(1/(2*(*tau))))*(((domain_constants.e[d][i]-current_node->u[d])*3)+(domain_constants.e[d][i]*9*((domain_constants.e[0][i]*current_node->u[0])+(domain_constants.e[1][i]*current_node->u[1]))));
		#endif
		}

		force_term[i] = 0;
		#pragma unroll
		for(d = 0; d<DIM; d++)
		{
			force_term[i] += F_coeff[d]*current_node->F[d];
		}

		eu = e_mul_u(current_node, &i);
		f_eq[i] = current_node->rho*domain_constants.omega[i]*(1.0+eu+(0.5*eu*eu)-u_sq);
	}

	if (current_node->c_smag>0) turbulent_viscosity(current_node, f_eq, tau);

	for(int i =0;i<Q;i++)
	{
		collision_bgk = (1.0/(*tau)) * (current_node->f[i]-f_eq[i]);
		collision_s = current_node->f[domain_constants.opp[i]]-current_node->f[i];

		tmp[i] = current_node->f[i] - (1-(current_node->B))*(collision_bgk) + (current_node->B)*collision_s + (1-(current_node->B))*force_term[i];
	}

	for(int i =0;i<Q;i++)
	{
		current_node->f[i] = tmp[i];
	}
}

__device__ void mrt_collision(Node *current_node, double *tau)
{
	double m_eq[Q],m[Q];
	
	#ifdef D2Q9
		meq_d2q9(current_node,m_eq);
	#endif

	#ifdef D3Q15
		meq_d3q15(current_node, m_eq);
	#endif

	for(int i = 0; i<Q; i++)
	{
		m[i] = 0;
		for(int j=0; j<Q; j++)
		{
			m[i] = m[i] + domain_constants.M[i][j]*current_node->f[j];
		}
	}

	for(int i = 0; i<Q;i++)
	{
		m[i] = domain_constants.tau_mrt[i]*(m[i]-m_eq[i]);
	}

	for(int i = 0; i<Q; i++)
	{
		//reuse m_eq to save on memory...
		m_eq[i] = 0;
		for(int j=0; j<Q; j++)
		{
			
			m_eq[i] = m_eq[i] + domain_constants.M_inv[i][j]*m[j];
		}
		current_node->f[i] = current_node->f[i] - m_eq[i]; // m_eq here is not equilibrium distribution, 
	}													   // it is the result of previous computation!!!!
	
}

__device__ void mrt_guo_collision(Node *current_node, double *tau)
{
	double m_eq[Q],m[Q], F_coeff[DIM], force_term[Q];
	int d;
	
	// Add force contribution to velocity
	#pragma unroll
	for(d = 0; d<DIM; d++)
	{
		current_node->u[d] = current_node->u[d] + (1/(2*current_node->rho))*current_node->F[d];
	}

	// Calculate forcing term
	for(int i=0;i<Q;i++)
	{
		#pragma unroll
		for(d = 0; d<DIM; d++)
		{
		#if DIM > 2
			F_coeff[d] = domain_constants.omega[i]*(1-(1/(2*(*tau))))*(((domain_constants.e[d][i]-current_node->u[d])*3)+(domain_constants.e[d][i]*9*((domain_constants.e[0][i]*current_node->u[0])+(domain_constants.e[1][i]*current_node->u[1])+(domain_constants.e[2][i]*current_node->u[2]))));
		#else
			F_coeff[d] = domain_constants.omega[i]*(1-(1/(2*(*tau))))*(((domain_constants.e[d][i]-current_node->u[d])*3)+(domain_constants.e[d][i]*9*((domain_constants.e[0][i]*current_node->u[0])+(domain_constants.e[1][i]*current_node->u[1]))));
		#endif
		}

		force_term[i] = 0;
		#pragma unroll
		for(d = 0; d<DIM; d++)
		{
			force_term[i] += F_coeff[d]*current_node->F[d];
		}
	}

	// Calculate equilibrium distribution
	#ifdef D2Q9
		meq_d2q9(current_node,m_eq);
	#endif

	#ifdef D3Q15
		meq_d3q15(current_node, m_eq);
	#endif

	// Convert PDF's to MDF's (Momentum Distribution Function)
	for(int i = 0; i<Q; i++)
	{
		m[i] = 0;
		for(int j=0; j<Q; j++)
		{
			m[i] = m[i] + domain_constants.M[i][j]*current_node->f[j];
		}
	}

	// Execute MRT collision
	for(int i = 0; i<Q;i++)
	{
		m[i] = domain_constants.tau_mrt[i]*(m[i]-m_eq[i]);
	}

	// convert MDF's back to PDF's adding the result of collision and forcing
	for(int i = 0; i<Q; i++)
	{
		//reuse m_eq to save on memory...
		m_eq[i] = 0;
		for(int j=0; j<Q; j++)
		{
			
			m_eq[i] = m_eq[i] + domain_constants.M_inv[i][j]*m[j];
		}
		current_node->f[i] = current_node->f[i] - m_eq[i] + force_term[i]; // m_eq here is not equilibrium distribution, 
	}																		// it is the result of previous computation!!!!
	
}

__device__ void mrt_ntpor_collision(Node *current_node, double *tau)
{
	double m_eq[Q],m[Q], collision_s;
	
	#ifdef D2Q9
		meq_d2q9(current_node,m_eq);
	#endif

	#ifdef D3Q15
		meq_d3q15(current_node, m_eq);
	#endif

	for(int i = 0; i<Q; i++)
	{
		m[i] = 0;
		for(int j=0; j<Q; j++)
		{
			m[i] = m[i] + domain_constants.M[i][j]*current_node->f[j];
		}
	}

	for(int i = 0; i<Q;i++)
	{
		m[i] = domain_constants.tau_mrt[i]*(m[i]-m_eq[i]);
	}

	for(int i = 0; i<Q; i++)
	{
		//reuse m_eq to save on memory...
		m_eq[i] = 0;
		for(int j=0; j<Q; j++)
		{
			m_eq[i] = m_eq[i] + domain_constants.M_inv[i][j]*m[j];
		}

		collision_s = current_node->f[domain_constants.opp[i]]-current_node->f[i];

		current_node->f[i] = current_node->f[i] - (1-(current_node->B))*m_eq[i] + (current_node->B)*collision_s; // m_eq here is not equilibrium distribution, 
	}													   // it is the result of previous computation!!!!
	
}

__device__ void mrt_ntpor_guo_collision(Node *current_node, double *tau)
{
	double m_eq[Q],m[Q], F_coeff[DIM], force_term[Q], collision_s;
	int d;
	
	// Add force contribution to velocity
	#pragma unroll
	for(d = 0; d<DIM; d++)
	{
		current_node->u[d] = current_node->u[d] + (1/(2*current_node->rho))*current_node->F[d];
	}

	// Calculate forcing term
	for(int i=0;i<Q;i++)
	{
		#pragma unroll
		for(d = 0; d<DIM; d++)
		{
		#if DIM > 2
			F_coeff[d] = domain_constants.omega[i]*(1-(1/(2*(*tau))))*(((domain_constants.e[d][i]-current_node->u[d])*3)+(domain_constants.e[d][i]*9*((domain_constants.e[0][i]*current_node->u[0])+(domain_constants.e[1][i]*current_node->u[1])+(domain_constants.e[2][i]*current_node->u[2]))));
		#else
			F_coeff[d] = domain_constants.omega[i]*(1-(1/(2*(*tau))))*(((domain_constants.e[d][i]-current_node->u[d])*3)+(domain_constants.e[d][i]*9*((domain_constants.e[0][i]*current_node->u[0])+(domain_constants.e[1][i]*current_node->u[1]))));
		#endif
		}

		force_term[i] = 0;
		#pragma unroll
		for(d = 0; d<DIM; d++)
		{
			force_term[i] += F_coeff[d]*current_node->F[d];
		}
	}

	// Calculate equilibrium distribution
	#ifdef D2Q9
		meq_d2q9(current_node,m_eq);
	#endif

	#ifdef D3Q15
		meq_d3q15(current_node, m_eq);
	#endif

	// Convert PDF's to MDF's (Momentum Distribution Function)
	for(int i = 0; i<Q; i++)
	{
		m[i] = 0;
		for(int j=0; j<Q; j++)
		{
			m[i] = m[i] + domain_constants.M[i][j]*current_node->f[j];
		}
	}

	// Execute MRT collision
	for(int i = 0; i<Q;i++)
	{
		m[i] = domain_constants.tau_mrt[i]*(m[i]-m_eq[i]);
	}

	// convert MDF's back to PDF's adding the result of collision and forcing
	for(int i = 0; i<Q; i++)
	{
		//reuse m_eq to save on memory...
		m_eq[i] = 0;
		for(int j=0; j<Q; j++)
		{
			
			m_eq[i] = m_eq[i] + domain_constants.M_inv[i][j]*m[j];
		}

		collision_s = current_node->f[domain_constants.opp[i]]-current_node->f[i];

		current_node->f[i] = current_node->f[i] - (1-(current_node->B))*m_eq[i] + (1-(current_node->B))*force_term[i] + (current_node->B)*collision_s; // m_eq here is not equilibrium distribution, 
	}																		// it is the result of previous computation!!!!
	
}

__device__ void bounceback(Node *current_node, double *tau)
{
	double tmp[Q];
	for(int i=0;i<Q;i++)
	{
		tmp[i] = current_node->f[i];
	}

	for(int i=0;i<Q;i++)
	{
		current_node->f[i] = tmp[domain_constants.opp[i]];
	}

	current_node->u[0] = 0;
	current_node->u[1] = 0;
	#if DIM > 2
		current_node->u[2] = 0;
	#endif

	current_node->rho = 0;
}

__device__ void turbulent_viscosity(Node *current_node, double *f_eq, double *tau)
{
	double q_bar[DIM][DIM];
	double q_hat = 0.;

	for(int i = 0; i<DIM; i++)
	{
		for(int j = 0; j<DIM; j++)
		{
			for(int q = 0; q<Q; q++)
			{
				q_bar[i][j] = q_bar[i][j]+((double)domain_constants.e[i][q]*(double)domain_constants.e[j][q]*(current_node->f[q]-f_eq[q]));
			}
			q_hat = q_hat + sqrt((double)2*q_bar[i][j]*q_bar[i][j]);
		}
	}
	
	*tau = *tau+0.5*(sqrt(((*tau)*(*tau))+(2*sqrt((double)2)*(current_node->c_smag*current_node->c_smag)*(1/(current_node->rho*POW4(1/sqrt((double)3))))*q_hat))-*tau);
}

__device__ void meq_d2q9(Node *current_node, double *meq)
{
	double jx = current_node->rho*current_node->u[0];
	double jy = current_node->rho*current_node->u[1];

	meq[0] = current_node->rho;
	meq[1] = (-2*current_node->rho)+3*(jx*jx+jy*jy);
	meq[2] = current_node->rho-3*(jx*jx+jy*jy);
	meq[3] = jx;
	meq[4] = -jx;
	meq[5] = jy;
	meq[6] = -jy;
	meq[7] = jx*jx-jy*jy;
	meq[8] = jx*jy;
}

__device__ void meq_d3q15(Node *current_node, double *meq)
{
	double jx = current_node->rho*current_node->u[0];
	double jy = current_node->rho*current_node->u[1];
	double jz = current_node->rho*current_node->u[2];

	meq[0] = current_node->rho;
	meq[1] = (-1*current_node->rho)+(jx*jx+jy*jy+jz*jz);
	//meq[2] = -current_node->rho;
	meq[2] = current_node->rho-5*(jx*jx+jy*jy+jz*jz);
	meq[3] = jx;
	meq[4] = (-7/3)*jx;
	meq[5] = jy;
	meq[6] = (-7/3)*jy;
	meq[7] = jz;
	meq[8] = (-7/3)*jz;
	meq[9] = (2*jx*jx-(jy*jy+jz*jz));
	meq[10] = jy*jy-jz*jz;
	meq[11] = jx*jy;
	meq[12] = jy*jz;
	meq[13] = jx*jz;
	meq[14] = 0;
}
#endif
