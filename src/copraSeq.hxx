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
 * Move each vertex to its best community.
 * @param vcs communities vertex u is linked to (updated)
 * @param vcout total edge weight from vertex u to community C (updated)
 * @param vcom community each vertex belongs to (updated)
 * @param x original graph
 * @param vdom community each vertex belonged to
 * @returns number of changed vertices
 */
template <bool STRICT=false, class G, class K, class V, class FA, class FP>
K copraMoveIteration(vector<K>& vcs, vector<V>& vcout, vector<K>& vcom, const G& x, FA fa, FP fp) {
  K a = K();
  x.forEachVertexKey([&](auto u) {
    if (!fa(u)) return;
    K d = vcom[u];
    copraClearScan(vcs, vcout);
    copraScanCommunities(vcs, vcout, x, u, vcom);
    auto [c, w] = copraChooseCommunity<STRICT>(x, u, vcom, vcs, vcout);
    if (c && c!=d) { vcom[u] = c; ++a; fp(u); }
  });
  return a;
}




// COPRA-SEQ
// ---------

template <bool STRICT=false, class G, class K, class FA, class FP>
CopraResult<K> copraSeq(const G& x, const vector<K>* q, const CopraOptions& o, FA fa, FP fp) {
  using V = typename G::edge_value_type;
  int l = 0;
  K S = x.span();
  K N = x.order();
  vector<K> vcom(S), vcs;
  vector<V> vcout(S);
  float t = measureDuration([&]() {
    copraInitialize(vcom, x);
    for (l=0; l<o.maxIterations;) {
      K n = copraMoveIteration<STRICT>(vcs, vcout, vcom, x, fa, fp); ++l;
      PRINTFD("copraSeq(): l=%d, n=%d, N=%d, n/N=%f\n", l, n, N, float(n)/N);
      if (float(n)/N <= o.tolerance) break;
    }
  }, o.repeat);
  return {vcom, l, t};
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
  const vector<K>& vcom = *q;
  auto vaff = copraAffectedVerticesDeltaScreening<STRICT>(x, deletions, insertions, vcom);
  auto fa   = [&](auto u) { return vaff[u]==true; };
  return copraSeq<STRICT>(x, q, o, fa);
}




// COPRA-SEQ-DYNAMIC-FRONTIER
// --------------------------

template <bool STRICT=false, class G, class K, class V>
inline CopraResult<K> copraSeqDynamicFrontier(const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>* q, const CopraOptions& o={}) {
  K S = x.span();
  const vector<K>& vcom = *q;
  auto vaff = copraAffectedVerticesFrontier(x, deletions, insertions, vcom);
  auto fa = [&](auto u) { return vaff[u]==true; };
  auto fp = [&](auto u) { x.forEachEdgeKey(u, [&](auto v) { vaff[v] = true; }); };
  return copraSeq<STRICT>(x, q, o, fa, fp);
}
