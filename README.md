Single-threaded CPU-based Community OVerlap PRopagation Algorithm ([COPRA]) for
[community detection].

This is an implementation of a label-propagation based community detection
algorithm called **Community OVerlap PRopagation Algorithm (COPRA)**. Unlike
**RAK**, this algorithm uses *multiple labels per vertex*, with each label having an
associated *belonging coefficient* (which sums to `1`). The algorithm is as follows.

1. Each vertex initializes as its own community (belonging=1).
2. Each iteration, a vertex collects labels from its neighborhood.
3. I excludes a vertex's own labels, although not explicitly mentioned in paper.
4. The collected labels are scaled by edge weights for weighted graph.
5. Each vertex picks labels above a certain threshold.
6. This threshold is inversely proportional to the max. number of labels.
7. If all labels are below threshold, pick a random best label.
8. I make a vertex join its own community if it has no labels (not mentioned).
9. Selected labels are normalized such that belonging coefficient sums to 1.
10. Repeat from 2 until convergence.

The authors of this paper mention to use a *mimimum vertices per community count*
to detect convergence. However i do not find it to be helpful. I instead use a
similar convergence condition as *RAK*, that is **count the number of vertices that**
**change thier best label**. Once this count falls below a certain fraction
(tolerance), i consider the algorithm to have converged. The authors also
mention using a *sychrounous* version of the algorithm, where labels of each
vertex is dependent only upon labels in previous iteration. However, i find
**asynchronous** approach to be converge faster (labels of each vertex can be
dependent upon labels in current iteration). As we focus of finding disjoint
communities, i consider the **best label of each vertex** as the final result.

[![](https://i.imgur.com/6UOli7q.png)][sheetp]

[![](https://i.imgur.com/7RUqa6l.png)][sheetp]

<br>
<br>

I vary the **tolerance** from `0.1` to `0.0001` (as with [RAK]), and adjust the **max.**
**number of labels** from `1` to `32`.

[![](https://i.imgur.com/70hxdGa.png)][sheetp]

[![](https://i.imgur.com/sG6ASLc.png)][sheetp]

[![](https://i.imgur.com/2NnJXny.png)][sheetp]

<br>
<br>

On average, it *seems that using a single label is best* in terms of time as well
as modularity. This is particulary true in case of road networks, but *not so in*
*case of other classes of graphs* (web graphs, social networks, collaboration
networks). For example, *web graphs* such as `web-Stanford` and `web-BerkStan` achieve
best modularity with **max. labels** of `8`, `web-Google` does best with **max. labels** of
`32`, and `web-NotreDame` does best with **max. labels** of `16`. **Max. labels** of `4`-`16`
would be a good choice for such graphs.

[![](https://i.imgur.com/zqq3eJO.png)][sheetp]

<br>
<br>

In addition it seems that on average, *making the tolerance tighter than 0.01 has*
*no beneficial effect on modularity*. However, *tighter tolerance does not help*
*with* graphs such as `web-NotreDame`, `coAuthorsDBLP`, and *social networks*. It seems
a **tolerance** of `0.01` would be a good choice on average.

[![](https://i.imgur.com/4rd819J.png)][sheetp]

<br>
<br>

Both **RAK** and **COPRA** approaches can obtain *disconnected communities*. This issue
can be resolved by splitting such communities into separate communities in a
*post-processing step*. I do **not** do that here. We may do it when comparing these
approaches.

All outputs are saved in a [gist] and a small part of the output is listed here.
Some [charts] are also included below, generated from [sheets]. The input data
used for this experiment is available from the [SuiteSparse Matrix Collection].
This experiment was done with guidance from [Prof. Kishore Kothapalli] and
[Prof. Dip Sankar Banerjee].


[COPRA]: https://arxiv.org/abs/0910.5516
[RAK]: https://arxiv.org/abs/0709.2938
[community detection]: https://en.wikipedia.org/wiki/Community_search
[previous experiment]: https://github.com/puzzlef/rak-communities-seq
[Prof. Dip Sankar Banerjee]: https://sites.google.com/site/dipsankarban/
[Prof. Kishore Kothapalli]: https://faculty.iiit.ac.in/~kkishore/
[SuiteSparse Matrix Collection]: https://sparse.tamu.edu

<br>

```bash
$ g++ -std=c++17 -O3 main.cxx
$ ./a.out ~/data/web-Stanford.mtx
$ ./a.out ~/data/web-BerkStan.mtx
$ ...

# Loading graph /home/subhajit/data/web-Stanford.mtx ...
# order: 281903 size: 2312497 [directed] {}
# order: 281903 size: 3985272 [directed] {} (symmetricize)
# [-0.000497 modularity] noop
# [00206.635 ms; 0003 iters.; 0.836576819 modularity] copraSeqStatic {labels=01, tolerance=1e-01}
# [00248.339 ms; 0003 iters.; 0.838712633 modularity] copraSeqStatic {labels=02, tolerance=1e-01}
# [00369.730 ms; 0003 iters.; 0.875163972 modularity] copraSeqStatic {labels=04, tolerance=1e-01}
# [00550.753 ms; 0003 iters.; 0.878541648 modularity] copraSeqStatic {labels=08, tolerance=1e-01}
# [00675.437 ms; 0003 iters.; 0.881162941 modularity] copraSeqStatic {labels=16, tolerance=1e-01}
# [01130.822 ms; 0003 iters.; 0.866800904 modularity] copraSeqStatic {labels=32, tolerance=1e-01}
# [00206.794 ms; 0003 iters.; 0.836576819 modularity] copraSeqStatic {labels=01, tolerance=5e-02}
# [00243.418 ms; 0003 iters.; 0.838712633 modularity] copraSeqStatic {labels=02, tolerance=5e-02}
# [00475.227 ms; 0004 iters.; 0.887068212 modularity] copraSeqStatic {labels=04, tolerance=5e-02}
# [00834.529 ms; 0005 iters.; 0.892077446 modularity] copraSeqStatic {labels=08, tolerance=5e-02}
# [01010.437 ms; 0005 iters.; 0.897701085 modularity] copraSeqStatic {labels=16, tolerance=5e-02}
# [01706.960 ms; 0005 iters.; 0.885300398 modularity] copraSeqStatic {labels=32, tolerance=5e-02}
# ...
```

<br>
<br>


## References

- [Finding overlapping communities in networks by label propagation; Steve Gregory (2010)](https://iopscience.iop.org/article/10.1088/1367-2630/12/10/103018)
- [Near linear time algorithm to detect community structures in large-scale networks; Usha Nandini Raghavan et al. (2007)](https://arxiv.org/abs/0709.2938)
- [The University of Florida Sparse Matrix Collection; Timothy A. Davis et al. (2011)](https://doi.org/10.1145/2049662.2049663)
- [How to import VSCode keybindings into Visual Studio?](https://stackoverflow.com/a/62417446/1413259)
- [Configure X11 Forwarding with PuTTY and Xming](https://www.centlinux.com/2019/01/configure-x11-forwarding-putty-xming-windows.html)
- [Installing snap on CentOS](https://snapcraft.io/docs/installing-snap-on-centos)

<br>
<br>


[![](https://i.imgur.com/7GLy9tb.jpg)](https://www.youtube.com/watch?v=L-ZBWLYGSuY)<br>
[![ORG](https://img.shields.io/badge/org-puzzlef-green?logo=Org)](https://puzzlef.github.io)


[gist]: https://gist.github.com/wolfram77/417cfff0ee5e5b233056283fc42f78ac
[charts]: https://imgur.com/a/zmjeCyw
[sheets]: https://docs.google.com/spreadsheets/d/1jNoY9zpiMpmFuRLMSTLY--XArqXmH0oQbguMtWv_j9s/edit?usp=sharing
[sheetp]: https://docs.google.com/spreadsheets/d/e/2PACX-1vQ9ym6LJBJZkwYEYk3sjAPJEt0QT67UZCxmnnOibZYYREwGXUmXL4LvurXAme5MvHlKMXaX-DOLX9Js/pubhtml
