#pragma once

#include "TopologicalFaceSet.hh"

namespace OpenVolumeMesh {
/* Computes the 'Link' of a vertex in the topological sense */
    template<class _MESH>
    TopologicalFaceSet link(const _MESH &mesh,
                            const OpenVolumeMesh::VertexHandle &vertex);


/* Computes the 'Link' of an edge in the topological sense */
    template<class _MESH>
    TopologicalFaceSet link(const _MESH &mesh,
                            const OpenVolumeMesh::EdgeHandle &edge);


/* Returns the set of topological faces that are in the intersection of
the edge's vertices' links but not in the edge's link */
    template<class _MESH>
    TopologicalFaceSet link_outsiders(const _MESH &mesh,
                                      const OpenVolumeMesh::EdgeHandle &edge);


/* Checks that the 'link condition' for a particular edge holds.
i.e. If edge goes from vertex a to vertex b, it checks whether

Lk(a) [intersection] Lk(b) = Lk(ab)

is true or not.

This can typically be used to check if this edge can be safely collapsed
without hurting the mesh's topology. */
    template<class _MESH>
    bool link_condition(const _MESH &mesh,
                        const OpenVolumeMesh::EdgeHandle &eh);


    template<class _MESH>
    bool link_condition(const _MESH &mesh,
                        const OpenVolumeMesh::HalfEdgeHandle &heh);
}

#include <TopologicalLinkT_impl.hh>
