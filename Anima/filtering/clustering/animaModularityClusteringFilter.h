#pragma once

#include <vector>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

namespace anima {
	
	template <class DataType>
    class ModularityClusteringFilter
    {
	public:
		typedef vnl_matrix <DataType> DataHolderType;
		typedef std::vector < unsigned int > MembershipType;
		typedef std::vector < MembershipType > MembershipVectorType;
		typedef vnl_vector < double > VectorType;
		typedef vnl_matrix < double > MatrixType;
		
		ModularityClusteringFilter();
        virtual ~ModularityClusteringFilter();
		
		void SetInputData(DataHolderType &data);
		void Update();
		
		unsigned int GetClassMembership(unsigned int i) {return m_ClassMemberships[i];}
		MembershipType& GetClassMemberships() {return m_ClassMemberships;}
		MembershipType& GetReverseClassMembership(unsigned int i) {return m_ReverseClassMemberships[i];}
		MembershipVectorType& GetReverseClassMemberships() {return m_ReverseClassMemberships;}
		
		unsigned int GetNumberOfClusters() {return m_NumberOfClusters;}
        
	private:
		void InitializeClassMemberships();
		void GetModularityMatrix();
		void GetSubModularityMatrix(unsigned int parentIndex);
		double SpectralClustering(unsigned int parentIndex);
		void SplitCluster(unsigned int parentIndex);
		
		MembershipType m_ClassMemberships, m_Membership, m_SubMembership1, m_SubMembership2;
		MembershipVectorType m_ReverseClassMemberships;
		MatrixType m_ConnectivityMatrix, m_ModularityMatrix, m_SubModularityMatrix;
		VectorType m_DegreeMatrix, m_SubDegreeMatrix;
		
		std::vector < bool > m_PreventClusterFromSplitting;
		
		unsigned int m_NbInputs;
		unsigned int m_NumberOfClusters;
		
		double m_TotalNumberOfEdges;
		
		static const double m_ZeroThreshold;
	};
	
} // end namespace anima

#include "animaModularityClusteringFilter.hxx"


