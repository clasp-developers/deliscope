(require :asdf)
(asdf:load-asd (pathname "/home/meister/Development/clasp-main/extensions/seqan/src/seqan.asd"))
(asdf:load-asd (pathname "/home/meister/Development/analyze-del/src/analyze-del.asd"))
(asdf:load-system :analyze-del :force t)

(ad:load-encoding "/home/meister/Development/analyze-talon/del-codes")

(progn
  (defparameter short-file
    (sa:make-seq-file-in "/home/meister/Development/analyze-talon/pass-s2-full.fastq"))
  (defparameter title (sa:make-string :char-string))
  (defparameter seq (sa:make-string :dna5q-string))
  (defparameter align (sa:make-align :dna5q-string))
  )

(progn
  (defparameter short-file
    (sa:make-seq-file-in "/home/meister/Development/analyze-talon/fail-s2-full.fastq"))
  (defparameter title (sa:make-string :char-string))
  (defparameter seq (sa:make-string :dna5q-string))
  (defparameter align (sa:make-align :dna5q-string))
  )

(progn
  (sa:read-record title seq short-file)
  (format t "sequence: ~a~%" (sa:to-string seq))
  (analyze-del::analyze-sequence align seq :verbose t)
  )


