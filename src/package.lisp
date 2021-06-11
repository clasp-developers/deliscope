
(defpackage :analyze-del
  (:use :common-lisp)
  (:nicknames :ad)
  (:shadow sequence)
  (:export
   #:*minimum-phred-quality*
   #:sequences
   #:set-filter
   #:library-codes
   #:filter-seq-file
   #:load-encoding
   #:analyze
   #:parallel-analyze
   #:serial-analyze
   #:clear-analysis
   #:remove-rare-sequences
   #:join-codons
   #:filter-sequences
   #:bead-specific-counts
   #:sort-hits
   #:*analysis*
   #:save-results
   #:load-results
   #:multiple-redundant
   #:merge-hits
   ))
   
