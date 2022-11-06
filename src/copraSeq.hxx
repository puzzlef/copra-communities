#pragma once
#include <utility>
#include <vector>
#include <algorithm>
#include "_main.hxx"
#include "vertices.hxx"
#include "edges.hxx"
#include "csr.hxx"
#include "copra.hxx"

using std::tuple;
using std::vector;
using std::swap;




// COPRA-MOVE-ITERATION
// --------------------

/**
 * Move each vertex to its best community set.
 * @param vcs communities vertex u is linked to (updated)
 * @param vcout total edge weight from vertex u to community C (updated)
 * @param vcom community set each vertex belongs to (updated)
 * @param x original graph
 * @param vdom community set each vertex belonged to
 * @param vtot total edge weight of each vertex
 * @param B belonging coefficient threshold
 * @returns number of changed vertices
 */
template <bool STRICT=false, class G, class K, class V, size_t L, class FA, class FP>
K copraMoveIteration(vector<K>& vcs, vector<V>& vcout, vector<Labelset<K, V, L>>& vcom, const G& x, const vector<Labelset<K, V, L>>& vdom, const vector<V>& vtot, V B, FA fa, FP fp) {
  K a = K();
  x.forEachVertexKey([&](auto u) {
    if (!fa(u)) return;
    K d = vdom[u][0].first;
    copraClearScan(vcs, vcout);
    copraScanCommunities(vcs, vcout, x, u, vdom);
    copraSortScan<STRICT>(vcs, vcout);
    copraChooseCommunity(vcom[u], u, vcs, vcout, B*vtot[u]);
    K c = vcom[u][0].first;
    if (c!=d) { ++a; fp(u); }
  });
  return a;
}




// COPRA-SEQ
// ---------

template <size_t LABELS=8, bool STRICT=false, bool ASYNC=false, class G, class K, class FA, class FP>
CopraResult<K> copraSeq(const G& x, const vector<K>* q, const CopraOptions& o, FA fa, FP fp) {
  using V = typename G::edge_value_type;
  const size_t L = LABELS;
  int l = 0;
  K S = x.span();
  K N = x.order();
  V B = V(1)/o.maxMembership;
  vector<K> vcs;
  vector<V> vcout(S), vtot(S);
  vector<Labelset<K, V, L>> vcom(S), vdom(S);
  float t = measureDuration([&]() {
    copraVertexWeights(vtot, x);
    copraInitialize(vdom, x);
    for (l=0; l<o.maxIterations;) {
      K n = copraMoveIteration<STRICT>(vcs, vcout, ASYNC? vdom : vcom, x, vdom, vtot, B, fa, fp); ++l;
      PRINTFD("copraSeq(): l=%d, n=%d, N=%d, n/N=%f\n", l, n, N, float(n)/N);
      if (float(n)/N <= o.tolerance) break;
      if (!ASYNC) swap(vdom, vcom);
    }
  }, o.repeat);
  return {copraBestCommunities(vcom), l, t};
}
template <bool STRICT=false, class G, class K, class FA>
inline CopraResult<K> copraSeq(const G& x, const vector<K>* q, const CopraOptions& o, FA fa) {
  auto fp = [](auto u) {};
  return copraSeq<STRICT>(x, q, o, fa, fp);
}
template <bool STRICT=false, class G, class K>
inline CopraResult<K> copraSeq(const G& x, const vector<K>* q, const CopraOptions& o) {
  auto fa = [](auto u) { return true; };
  return copraSeq<STRICT>(x, q, o, fa);
}




// COPRA-SEQ-STATIC
// ----------------

template <bool STRICT=false, class G, class K>
inline CopraResult<K> copraSeqStatic(const G& x, const vector<K>* q=nullptr, const CopraOptions& o={}) {
  return copraSeq<STRICT>(x, q, o);
}




// COPRA-SEQ-DYNAMIC-DELTA-SCREENING
// ---------------------------------

template <bool STRICT=false, class G, class K, class V>
inline CopraResult<K> copraSeqDynamicDeltaScreening(const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>* q, const CopraOptions& o={}) {
  K S = x.span();
  const vector<Labelset<K, V, L>>& vcom = *q;
  auto vaff = copraAffectedVerticesDeltaScreening<STRICT>(x, deletions, insertions, vcom);
  auto fa   = [&](auto u) { return vaff[u]==true; };
  return copraSeq<STRICT>(x, q, o, fa);
}




// COPRA-SEQ-DYNAMIC-FRONTIER
// --------------------------

template <bool STRICT=false, class G, class K, class V>
inline CopraResult<K> copraSeqDynamicFrontier(const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>* q, const CopraOptions& o={}) {
  K S = x.span();
  const vector<Labelset<K, V, L>>& vcom = *q;
  auto vaff = copraAffectedVerticesFrontier(x, deletions, insertions, vcom);
  auto fa = [&](auto u) { return vaff[u]==true; };
  auto fp = [&](auto u) { x.forEachEdgeKey(u, [&](auto v) { vaff[v] = true; }); };
  return copraSeq<STRICT>(x, q, o, fa, fp);
}