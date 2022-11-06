#pragma once
#include <limits>
#include <utility>
#include <vector>
#include "_main.hxx"

using std::pair;
using std::tuple;
using std::vector;
using std::numeric_limits;
using std::make_pair;
using std::move;
using std::get;




// COPRA-OPTIONS
// -------------

// Maximum community memberships per vertex.
#define COPRA_MAX_MEMBERSHIP 8


struct CopraOptions {
  int   repeat;
  float tolerance;
  int   maxIterations;

  CopraOptions(int repeat=1, float tolerance=0.05, int maxIterations=20) :
  repeat(repeat), tolerance(tolerance), maxIterations(maxIterations) {}
};




// COPRA-RESULT
// ------------

template <class K>
struct CopraResult {
  vector<K> membership;
  int   iterations;
  float time;

  CopraResult(vector<K>&& membership, int iterations=0, float time=0) :
  membership(membership), iterations(iterations), time(time) {}

  CopraResult(vector<K>& membership, int iterations=0, float time=0) :
  membership(move(membership)), iterations(iterations), time(time) {}
};




// LABELSET
// --------

template <class K, class V, size_t L>
using Labelset = array<pair<K, V>, L>;




// COPRA-INITIALIZE
// ----------------

/**
 * Find the total edge weight of each vertex.
 * @param vtot total edge weight of each vertex (updated)
 * @param x original graph
 */
template <class G, class V>
void copraVertexWeights(vector<V>& vtot, const G& x) {
  x.forEachVertexKey([&](auto u) {
    vtot[u] = V();
    x.forEachEdge(u, [&](auto v, auto w) { vtot[u] += w; });
  });
}


/**
 * Initialize communities such that each vertex is its own community.
 * @param vcom community set each vertex belongs to (updated)
 * @param x original graph
 */
template <class G, class K, class V, size_t L>
inline void copraInitialize(vector<Labelset<K, V, L>>& vcom, const G& x) {
  x.forEachVertexKey([&](auto u) { vcom[u] = {make_pair(u, V(1))}; });
}




// COPRA-CHOOSE-COMMUNITY
// ----------------------

/**
 * Scan an edge community connected to a vertex.
 * @param vcs communities vertex u is linked to (updated)
 * @param vcout total edge weight from vertex u to community C (updated)
 * @param u given vertex
 * @param v outgoing edge vertex
 * @param w outgoing edge weight
 * @param vcom community set each vertex belongs to
 */
template <bool SELF=false, class K, class V, size_t L>
inline void copraScanCommunity(vector<K>& vcs, vector<V>& vcout, K u, K v, V w, const vector<Labelset<K, V, L>>& vcom) {
  if (!SELF && u==v) return;
  for (const auto& [c, b] : vcom[v]) {
    if (!b) break;  // TODO? b -> c
    if (!vcout[c]) vcs.push_back(c);
    vcout[c] += w*b;
  }
}


/**
 * Scan communities connected to a vertex.
 * @param vcs communities vertex u is linked to (updated)
 * @param vcout total edge weight from vertex u to community C (updated)
 * @param x original graph
 * @param u given vertex
 * @param vcom community set each vertex belongs to
 */
template <bool SELF=false, class G, class K, class V, size_t L>
inline void copraScanCommunities(vector<K>& vcs, vector<V>& vcout, const G& x, K u, const vector<Labelset<K, V, L>>& vcom) {
  x.forEachEdge(u, [&](auto v, auto w) { copraScanCommunity<SELF>(vcs, vcout, u, v, w, vcom); });
}


/**
 * Sort communities scan data by total edge weight / belongingness.
 * @param vcs communities vertex u is linked to (updated)
 * @param vcout total edge weight from vertex u to community C (updated)
 */
template <class K, class V>
inline void copraSortScan(vector<K>& vcs, const vector<V>& vcout) {
  auto fl = [&](auto c, auto d) { return vcout[c] > vcout[d]; };
  sortValues(vcs, fl);
}


/**
 * Clear communities scan data.
 * @param vcs communities vertex u is linked to (updated)
 * @param vcout total edge weight from vertex u to community C (updated)
 */
template <class K, class V>
inline void copraClearScan(vector<K>& vcs, vector<V>& vcout) {
  for (K c : vcs)
    vcout[c] = V();
  vcs.clear();
}


/**
 * Choose connected community with most weight.
 * @param x original graph
 * @param u given vertex
 * @param vcom community each vertex belongs to
 * @param vcs communities vertex u is linked to
 * @param vcout total edge weight from vertex u to community C
 * @param W edge weight threshold above which communities are chosen
 * @returns [best community, best edge weight to community]
 */
template <class G, class K, class V, size_t L>
inline Labelset<K, V, L> copraChooseCommunity(const G& x, K u, const vector<Labelset<K, V, L>>& vcom, const vector<K>& vcs, const vector<V>& vcout, V W) {
  K n = K(); V w = V();
  Labelset<K, V, L> labs;
  // 1. Find labels above threshold, or best below threshold.
  for (K c : vcs) {
    if (n>K() && vcout[c]<W) break;
    labs[n++] = {c, vcout[c]};
    w += vcout[c];
  }
  // 2. Normalize labels, such that belonging coefficient sums to 1.
  for (K i=0; i<n; ++i)
    labs[i].second /= w;
  if (!n) labs[0] = make_pair(u, V(1));
  return labs;
}




// COPRA-COUNT-COMMUNITIES
// -----------------------

/**
 * Count number of vertices belonging to each community.
 * @param gcs list of communities (updated)
 * @param gcnum number of vertices belonging to respective community (updated)
 * @param x original graph
 * @param vcom community set each vertex belongs to
 */
template <class G, class K, class V, size_t L>
inline void copraCountCommunities(vector<K>& gcs, vector<K>& gcnum, const G& x, const vector<Labelset<K, V, L>>& vcom) {
  x.forEachVertexKey([&](auto u) {
    for (const auto& [c, b] : vcom[u]) {
      if (!b) break;
      if (!gcnum[c]) gcs.push_back(c);
      ++gcnum[c];
    }
  });
}


/**
 * Clear communities count data.
 * @param gcs list of communities (updated)
 * @param gcnum number of vertices belonging to respective community (updated)
 */
template <class K>
inline void copraClearCount(vector<K>& gcs, vector<K>& gcnum) {
  for (K c : gcs)
    gcnum[c] = K();
  gcs.clear();
}


/**
 * Get minimum community count.
 * @param gcs list of communities (updated)
 * @param gcnum number of vertices belonging to respective community (updated)
 */
template <class K>
inline K copraMinCount(vector<K> gcs, vector<K>& gcnum) {
  K min = numeric_limits<K>::max();
  for (K c : gcs)
    min = min <= gcnum[c]? min : gcnum[c];
  return min;
}




// COPRA-BEST-COMMUNITIES
// ----------------------

template <class K, class V, size_t L>
inline vector<K> copraBestCommunities(const vector<Labelset<K, V, L>>& vcom) {
  K S = vcom.size();
  vector<K> a(S);
  for (size_t i=0; i<S; ++i)
    a[i] = vcom[i][0].first;
  return a;
}




// COPRA-AFFECTED-VERTICES-DELTA-SCREENING
// ---------------------------------------
// Using delta-screening approach.
// - All edge batches are undirected, and sorted by source vertex-id.
// - For edge additions across communities with source vertex `i` and highest modularity changing edge vertex `j*`,
//   `i`'s neighbors and `j*`'s community is marked as affected.
// - For edge deletions within the same community `i` and `j`,
//   `i`'s neighbors and `j`'s community is marked as affected.

/**
 * Find the vertices which should be processed upon a batch of edge insertions and deletions.
 * @param x original graph
 * @param deletions edge deletions for this batch update (undirected, sorted by source vertex id)
 * @param insertions edge insertions for this batch update (undirected, sorted by source vertex id)
 * @param vcom community set each vertex belongs to
 * @param vtot total edge weight of each vertex
 * @param B belonging coefficient threshold
 * @returns flags for each vertex marking whether it is affected
 */
template <class FLAG=bool, class G, class K, class V, size_t L>
auto copraAffectedVerticesDeltaScreening(const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<Labelset<K, V, L>>& vcom, const vector<V>& vtot, V B) {
  K S = x.span();
  vector<K> vcs; vector<V> vcout(S);
  vector<FLAG> vertices(S), neighbors(S), communities(S);
  for (const auto& [u, v] : deletions) {
    K cu = vcom[u][0].first;
    K cv = vcom[v][0].first;
    if (cu!=cv) continue;
    vertices[u]  = true;
    neighbors[u] = true;
    communities[cv] = true;
  }
  for (size_t i=0; i<insertions.size();) {
    K u = get<0>(insertions[i]);
    copraClearScan(vcs, vcout);
    for (; i<insertions.size() && get<0>(insertions[i])==u; ++i) {
      K v  = get<1>(insertions[i]);
      V w  = get<2>(insertions[i]);
      K cu = vcom[u][0].first;
      K cv = vcom[v][0].first;
      if (cu==cv) continue;
      copraScanCommunity(vcs, vcout, u, v, w, vcom);
    }
    auto labs = copraChooseCommunity(x, u, vcom, vcs, vcout, B*vtot[u]);
    K cu = vcom[u][0].first;
    K cl = labs[0].first;
    if (cl==cu) continue;
    vertices[u]  = true;
    neighbors[u] = true;
    communities[cl] = true;
  }
  x.forEachVertexKey([&](auto u) {
    K cu = vcom[u][0].first;
    if (neighbors[u]) x.forEachEdgeKey(u, [&](auto v) { vertices[v] = true; });
    if (communities[cu]) vertices[u] = true;
  });
  return vertices;
}




// COPRA-AFFECTED-VERTICES-FRONTIER
// --------------------------------
// Using frontier based approach.
// - All source and destination vertices are marked as affected for insertions and deletions.
// - For edge additions across communities with source vertex `i` and destination vertex `j`,
//   `i` is marked as affected.
// - For edge deletions within the same community `i` and `j`,
//   `i` is marked as affected.
// - Vertices whose communities change in local-moving phase have their neighbors marked as affected.

/**
 * Find the vertices which should be processed upon a batch of edge insertions and deletions.
 * @param x original graph
 * @param deletions edge deletions for this batch update (undirected, sorted by source vertex id)
 * @param insertions edge insertions for this batch update (undirected, sorted by source vertex id)
 * @param vcom community set each vertex belongs to
 * @returns flags for each vertex marking whether it is affected
 */
template <class FLAG=bool, class G, class K, class V, size_t L>
auto copraAffectedVerticesFrontier(const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<Labelset<K, V, L>>& vcom) {
  K S = x.span();
  vector<FLAG> vertices(S);
  for (const auto& [u, v] : deletions) {
    K cu = vcom[u][0].first;
    K cv = vcom[v][0].first;
    if (cu!=cv) continue;
    vertices[u] = true;
  }
  for (const auto& [u, v, w] : insertions) {
    K cu = vcom[u][0].first;
    K cv = vcom[v][0].first;
    if (cu==cv) continue;
    vertices[u] = true;
  }
  return vertices;
}
