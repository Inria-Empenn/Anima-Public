#pragma once

#include <vector>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

namespace anima {
	
	template <class DataType, unsigned int DataDimension>
    class SPANClusteringFilter
    {
	public:
		typedef std::vector < DataType > DataHolderType;
		typedef std::vector < unsigned int > MembershipType;
		typedef std::vector < MembershipType > MembershipVectorType;
		typedef vnl_vector < double > VectorType;
		typedef vnl_matrix < double > MatrixType;
		
		SPANClusteringFilter();
        virtual ~SPANClusteringFilter();
		
		void SetInputData(DataHolderType &data);
		void SetOrientationData(bool orientationData) {m_OrientationData = orientationData;}
		void Update();
		
		unsigned int GetClassMembership(unsigned int i) {return m_ClassMemberships[i];}
		MembershipType& GetClassMemberships() {return m_ClassMemberships;}
		MembershipType& GetReverseClassMembership(unsigned int i) {return m_ReverseClassMemberships[i];}
		MembershipVectorType& GetReverseClassMemberships() {return m_ReverseClassMemberships;}
		
		unsigned int GetNumberOfClusters() {return m_NumberOfClusters;}
        
	private:
		void InitializeClassMemberships();
		MatrixType ComputeCMatrix(unsigned int parentClusterIndex, VectorType &BtPerOne);
		void FastBipartioning(MatrixType &data, unsigned int parentIndex, MembershipType &members1, MembershipType &members2);
		double GetNewmanGirvanModularity(MembershipType &members1, MembershipType &members2, VectorType &BtPerOne);
		void SplitCluster(unsigned int parentCluster);
		void GetConnectivityMatrix();
		void GetModularityMatrix();
		MatrixType GetSubModularityMatrix(unsigned int parentIndex);
		double SpectralClustering(unsigned int parentIndex, MembershipType &members1, MembershipType &members2);
		
		MembershipType m_ClassMemberships;
		MembershipVectorType m_ReverseClassMemberships;
		DataHolderType m_InputData;
		MatrixType m_ConnectivityMatrix, m_ModularityMatrix;
		VectorType m_DegreeMatrix;
		
		bool m_OrientationData;
		std::vector < bool > m_PreventClusterFromSplitting;
		
		unsigned int m_NbInputs;
		unsigned int m_NumberOfClusters;
		
		double m_TotalNumberOfEdges;
	};
	
} // end namespace anima

#include "animaSPANClusteringFilter.hxx"


