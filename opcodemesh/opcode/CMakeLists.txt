cmake_minimum_required(VERSION 2.6)

INCLUDE_DIRECTORIES (${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES ("${CMAKE_CURRENT_SOURCE_DIR}/Ice")

add_library(opcode STATIC
		   OPC_AABBCollider.cpp 
		   OPC_AABBTree.cpp 
		   OPC_BaseModel.cpp
		   OPC_BoxPruning.cpp 
		   OPC_Collider.cpp 
		   OPC_Common.cpp
		   OPC_HybridModel.cpp
		   OPC_LSSCollider.cpp 
		   OPC_MeshInterface.cpp
		   OPC_Model.cpp	
		   OPC_OBBCollider.cpp 
		   Opcode.cpp
		   OPC_OptimizedTree.cpp 
		   OPC_Picking.cpp 
		   OPC_PlanesCollider.cpp
		   OPC_RayCollider.cpp 
		   OPC_SphereCollider.cpp	
		   OPC_SweepAndPrune.cpp
		   OPC_TreeBuilders.cpp 
		   OPC_TreeCollider.cpp 
		   OPC_VolumeCollider.cpp
		   StdAfx.cpp 
		   Ice/IceAABB.cpp 
		   Ice/IceContainer.cpp
		   Ice/IceHPoint.cpp 
		   Ice/IceIndexedTriangle.cpp
		   Ice/IceMatrix3x3.cpp 
		   Ice/IceMatrix4x4.cpp 
		   Ice/IceOBB.cpp
		   Ice/IcePlane.cpp 
		   Ice/IcePoint.cpp 
		   Ice/IceRandom.cpp
		   Ice/IceRay.cpp 
		   Ice/IceRevisitedRadix.cpp 
		   Ice/IceSegment.cpp
		   Ice/IceTriangle.cpp 
		   Ice/IceUtils.cpp)

