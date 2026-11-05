// Research prototype ONLY: no production use, no robust I/O error semantics.
// Rcpp::sourceCpp with libBigWig includes and shared library supplied externally.
#include <Rcpp.h>
#include <bigWig.h>
#include <vector>
#include <algorithm>
#include <memory>
#include <cstdint>
using namespace Rcpp;
struct Exon { uint32_t lo,hi; R_xlen_t off; };
struct Span { uint32_t lo,hi; };
struct BwCloser { void operator()(bigWigFile_t* p) const {if(p) bwClose(p);} };
struct IvCloser { void operator()(bwOverlappingIntervals_t* p) const {if(p) bwDestroyOverlappingIntervals(p);} };
// [[Rcpp::export]]
NumericVector libbigwig_exons(std::string file, std::string chrom, IntegerVector starts,
                             IntegerVector ends, bool reverse_each, int gap=-1) {
  if(starts.size()!=ends.size() || gap < -1) stop("invalid arguments");
  static const int init_status = bwInit(1<<17);
  if (init_status) stop("libBigWig initialization failed");
  std::unique_ptr<bigWigFile_t, BwCloser> bw(bwOpen(file.c_str(),nullptr,"r"));
  if(!bw || bwGetTid(bw.get(),chrom.c_str())==UINT32_MAX) stop("cannot open file or chromosome");
  std::vector<Exon> exons; R_xlen_t total=0;
  for (R_xlen_t i=0;i<starts.size();i++) {
    if(starts[i]==NA_INTEGER || ends[i]==NA_INTEGER || starts[i]<1 || ends[i]<starts[i]) stop("invalid range");
    exons.push_back({static_cast<uint32_t>(starts[i]-1),static_cast<uint32_t>(ends[i]),total});
    total += static_cast<R_xlen_t>(ends[i])-starts[i]+1;
  }
  NumericVector out(total); std::vector<Span> spans;
  for(const auto& e:exons) spans.push_back({e.lo,e.hi});
  if(gap>=0) {
    std::sort(spans.begin(),spans.end(),[](const Span&a,const Span&b){return a.lo < b.lo;});
    std::vector<Span> merged;
    for(const auto& s:spans) {
      if(merged.empty() || static_cast<uint64_t>(s.lo) > static_cast<uint64_t>(merged.back().hi)+gap) merged.push_back(s);
      else merged.back().hi=std::max(merged.back().hi,s.hi);
    }
    spans.swap(merged);
  }
  for(const auto& span:spans) {
    std::unique_ptr<bwOverlappingIntervals_t, IvCloser> iv(bwGetOverlappingIntervals(bw.get(),chrom.c_str(),span.lo,span.hi));
    if(!iv) continue; // prototype cannot distinguish an empty interval from an I/O error
    for (const auto& e:exons) {
      // Independent exon queries use their matching exon only, including duplicates.
      if (e.lo<span.lo || e.hi>span.hi) continue;
      for(uint32_t j=0;j<iv->l;j++) {
        const uint32_t lo=std::max(e.lo,iv->start[j]), hi=std::min(e.hi,iv->end[j]);
        if(lo>=hi) continue;
        for(uint32_t p=lo;p<hi;p++) out[e.off+(reverse_each ? e.hi-1-p : p-e.lo)] = iv->value[j];
      }
    }
  }
  return out;
}
