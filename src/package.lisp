
(defpackage :deliscope
  (:use :common-lisp)
  (:nicknames :ds)
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



   #:compare
   #:create-sequence-parser
   #:filtered-library-codes
   #:summarize
   #:redundancies
   #:save-csv
   ))
   
