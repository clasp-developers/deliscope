(in-package :analyze-del)

(defparameter *simple-score* (sa:make-simple-score 0 -1 -2))
(defparameter *align-config* (sa:make-align-config :t-nil-nil-t))
(defparameter *del-code-format* "~d~a~2,'0d")

(defparameter *minimum-phred-quality* 20)
(defparameter *low-quality-max* 2)
(defparameter *minimum-counts* 100)

(defclass analysis ()
  ((file-names :initarg :file-names :accessor file-names)
   (sequences :initform (make-hash-table :test 'equal) :initarg :sequences :accessor sequences)
   (total-sequence-count :initform 0 :initarg :total-sequence-count :accessor total-sequence-count)
   (good-sequence-count :initform 0 :initarg :good-sequence-count :accessor good-sequence-count)))

(defmethod summarize (analysis)
  (format t "Files:~%~s~%" (file-names analysis))
  (format t "Total sequences: ~a~%" (total-sequence-count analysis))
  (format t "Good sequences:  ~a~%" (good-sequence-count analysis))
  (format t "Sample sequences:~%")
  (show-hash-table (sequences analysis)))

(defclass all-codes ()
  ((column-codes :initarg :column-codes :accessor column-codes)
   (forward-primer :initarg :forward-primer :accessor forward-primer)
   (end :initarg :end :accessor end)))

(defclass column-codes ()
  ((codes :initarg :codes :accessor codes)
   (start :initarg :start :accessor start)
   (end :initarg :end :accessor end)))

(defclass one-code ()
  ((del-name :initarg :del-name :accessor del-name)
   (del-index :initarg :del-index :accessor del-index)
   (sequence :initarg :sequence :accessor sequence)
   (code :initarg :code :accessor code)))

(defmethod print-object ((object one-code) stream)
  (print-unreadable-object (object stream :type t)
    (format stream "~a ~a ~a" (del-name object) (del-name object) (sa:to-string (code object)))))

(defclass match ()
  ((score :initarg :score :accessor score)
   (low-quality-count :initarg :low-quality-count :accessor low-quality-count)
   (quality-string :initarg :quality-string :accessor quality-string)
   (row-code :initarg :row-code :accessor row-code)))

(defmethod print-object ((object match) stream)
  (print-unreadable-object (object stream :type t)
    (format stream "~a ~a ~a ~s"
            (score object)
            (let ((*print-base* 16))
              (format nil "~a" (del-name (row-code object))))
            (low-quality-count object)
            (quality-string object))))


(defun hamming-distance (aaa bbb)
  (let ((hamming 0))
    (loop for pos below (length aaa)
          for ca = (elt aaa pos)
          for cb = (elt bbb pos)
          when (char/= ca cb)
            do (incf hamming)
          finally (return hamming))))


(defun hamming-matrix (seqs)
  (let* ((seq-len (length (first seqs)))
         (mat (make-array (list seq-len seq-len) :element-type 'fixnum)))
    (loop for posa from 0 below (1- seq-len)
          for seqa = (elt seqs posa)
          do (loop for posb from (1+ posa) below seq-len
                 for seqb = (elt seqs posb)
                   for hamming = (hamming-distance seqa seqb)
                   do (setf (aref mat posa posb) hamming)))
  mat))


(defun calculate-del-name (digit1 digit2 digit34 &key debug)
  (declare (ignore debug))
  (let* ((hex34 (if (= digit34 10) #x10 digit34))
         (res (+ (* digit1 #x1000) (* digit2 #x0100) hex34)))
      (let ((*print-base* 16))
        (format nil "~a" res))))

(defun join-del-codes (code1 code2)
  (concatenate 'string code1 code2))

(defun build-codes-for-column (left-overhang input-codes right-overhang pair-index sub-pair-index)
  (let* ((left-pos (+ (* pair-index 22) (* sub-pair-index 11)))
         (right-pos (+ 14 (* pair-index 22) (* sub-pair-index 11)))
         (codes (loop for index from 0 below (length input-codes)
                      for code in input-codes
                      for overhang-code = (concatenate 'string left-overhang code right-overhang)
                      for digit1 = (1+ sub-pair-index)
                      for digit2 = (+ 1 (* pair-index 2) sub-pair-index)
                      for digit34 = (1+ index)
                      collect (make-instance 'one-code
                                             :code (sa:make-string :dna5q-string overhang-code)
                                             :sequence overhang-code
                                             :del-name (calculate-del-name digit1 digit2 digit34)
                                             :del-index (1- digit34)))))
    (format t "~a ~a~%" left-pos right-pos)
    (make-instance 'column-codes
                   :codes codes
                   :start left-pos
                   :end right-pos)))

(defvar *all-codes* nil)

(defun load-encoding (del-code-file)
  (let* ((fin (open del-code-file :direction :input))
         forward-primer overhangs r1xseqs r2xseqs)
    (loop named reader
          for cmd = (read fin nil :eof)
          do (cond
               ((eq cmd :forward-primer) (setf forward-primer (read fin)))
               ((eq cmd :overhangs) (setf overhangs (read fin)))
               ((eq cmd :1xyy) (setf r1xseqs (read fin)))
               ((eq cmd :2xyy) (setf r2xseqs (read fin)))
               ((eq cmd :eof) (return-from reader nil))))
    (close fin)
    (let ((*print-pretty* t))
      (format t "1x sequences hamming distances~%~a~%" (hamming-matrix r1xseqs))
      (format t "2x sequences hamming distances~%~a~%" (hamming-matrix r2xseqs)))
    (format t "Codon offsets~%")
    (let ((column-codes (loop for cur = overhangs then (cdr cur)
                              for pair-index in '(0 0 1 1 2 2 3 3 4 4 )
                              for sub-pair-index in '(0 1 0 1 0 1 0 1 0 1)
                              for left-overhang = (car cur)
                              for right-overhang = (cadr cur)
                              for sequences in (list r1xseqs r2xseqs
                                                     r1xseqs r2xseqs
                                                     r1xseqs r2xseqs
                                                     r1xseqs r2xseqs
                                                     r1xseqs r2xseqs )
                              collect (build-codes-for-column left-overhang
                                                              sequences
                                                              right-overhang
                                                              pair-index
                                                              sub-pair-index))))
      (setf *all-codes* (make-instance 'all-codes
                                       :forward-primer (sa:make-string :dna5q-string forward-primer)
                                       :column-codes column-codes
                                       :end (end (car (last column-codes))))))))


(defparameter *library-code* nil)

(defun set-filter (&key library-code (minimum-counts 100))
  "Set the library code to use as a filter for sequences"
  (setf *library-code* library-code)
  (setf *minimum-counts* minimum-counts))

(defun analyze-column-using-align (align sequence column-codes start adjust low-quality-count quality-string &key verbose)
  (let* ((abs-start (+ start (start column-codes) adjust))
         (unsorted (loop for row-code in (codes column-codes)
                         for name = (del-name row-code)
                         for code = (code row-code)
                         for lower-diag = (- abs-start 1)
                         for upper-diag = (+ lower-diag 1) ;;+ start end-offset )
                         for score-code = (let (#+(or)(align (sa:make-align :dna5q-string)))
                                            (sa:resize (sa:rows& align) 2)
                                            (sa:assign-source (sa:row& align 0) sequence)
                                            (sa:assign-source (sa:row& align 1) code)
                                            (let ((score (sa:global-alignment align *simple-score* *align-config*
                                                                              lower-diag
                                                                              upper-diag )))
                                              (when verbose
                                                (format t "Score: ~a~%" score)
                                                (format t "~a~%" (sa:to-string align)))
                                              (make-instance 'match :score score
                                                                    :low-quality-count low-quality-count
                                                                    :quality-string quality-string
                                                                    :row-code row-code)))
                         collect score-code into uns
                         if (zerop (score score-code))
                           do (return uns)
                         finally (return uns))))
    (when verbose (format t "unsorted alignments: ~s~%" unsorted))
    (first (sort unsorted #'> :key #'score))))

(defun analyze-column-fast (sequence column-codes start adjust &key verbose)
  (let* ((seq-string (sa:to-string sequence))
         (abs-start (+ start (start column-codes) adjust))
         (abs-end (+ start (end column-codes) adjust)))
    (when verbose (format t "abs-start: ~d  abs-end ~d Subsequence: ~a~%" abs-start abs-end (subseq (sa:to-string sequence) (1- abs-start) (1- abs-end))))
    (loop for row-code in (codes column-codes)
          for code-seq = (sequence row-code)
          when (string= seq-string code-seq :start1 (1- abs-start) :end1 (1- abs-end))
            do (return-from analyze-column-fast (values t row-code)))
    (values nil)))

(defun analyze-column (align sequence column-codes start adjust &key verbose)
  (let* ((abs-start (+ start (start column-codes) adjust))
         (abs-end (+ start (end column-codes) adjust)))
    (multiple-value-bind (low-quality-count qual-string)
        (sa:count-quality-value-less-than sequence *minimum-phred-quality* abs-start abs-end)
      (multiple-value-bind (fast-worked row-code)
          (analyze-column-fast sequence column-codes start adjust :verbose verbose)
        (when fast-worked 
          (return-from analyze-column (make-instance 'match :score 0
                                                            :low-quality-count low-quality-count
                                                            :row-code row-code
                                                            :quality-string qual-string))))
      (analyze-column-using-align align sequence column-codes start adjust low-quality-count qual-string :verbose verbose))))

(defun analyze-forward-primer (align sequence &key verbose)
  (sa:resize (sa:rows& align) 2)
  (sa:assign-source (sa:row& align 0) sequence)
  (sa:assign-source (sa:row& align 1) (forward-primer *all-codes*))
  (let ((fwd-score (sa:global-alignment align *simple-score* *align-config*))
        (start (1+ (sa:to-view-position (sa:row& align 1) (sa:length (forward-primer *all-codes*))))))
    (when verbose (format t "Alignment: ~%~a~%" (sa:to-string align)))
    (values start fwd-score)))

(defun analyze-sequence (align sequence &key verbose)
  (multiple-value-bind (start fwd-score)
      (analyze-forward-primer align sequence :verbose verbose)
    (declare (ignore fwd-score))
    (when verbose
      (format t "len fwd: ~a~%" (sa:length (forward-primer *all-codes*)))
      (format t "3' end: ~a~%" start)
      (format t "~a~%" (sa:to-string align)))
    (when (> (sa:length sequence) (+ start (end *all-codes*)))
      (let ((result (loop for column-codes in (column-codes *all-codes*)
                          for start-offset = (start column-codes)
                          for end-offset = (end column-codes)
                          for best-digit = (analyze-column align sequence column-codes start 0 :verbose verbose)
                          do (when verbose (format t "offset: ~a best-digits: ~a~%" start-offset best-digit))
                          collect best-digit)))
        (when verbose
          (format t "result: ~a~%" result))
        result))))


(defvar *analysis* (make-hash-table :test 'equal))
(defvar *original-analysis* (make-instance 'analysis))
(defvar *common-sequences* (make-hash-table :test 'equal))
(defvar *bead-specific (make-hash-table :test 'equal))


(defun clear-analysis (&optional file-names)
  (setf *original-analysis* (make-instance 'analysis
                                           :file-names file-names))
  (clrhash *common-sequences*))

(defun accept-result (result)
  (and result
       (every (lambda (x) (>= (score x) -2)) result)
       (every (lambda (x) (< (low-quality-count x) *low-quality-max*)) result)))

  
(defun analyze-seq-file (file-name seq-file &key (num 10) verbose pass-file fail-file)
  (let ((title (sa:make-string :char-string))
        (seq (sa:make-string :dna5q-string))
        (align (sa:make-align :dna5q-string)))
    (loop named read-loop
          for index from 0
          for res = (progn
                      (when (sa:at-end seq-file) (return-from read-loop nil))
                      (sa:read-record title seq seq-file)
                      (analyze-sequence align seq :verbose verbose))
          unless (and num (< num 100000))
            do (when (= 0 (mod index 100000))
                 (format t "Done ~a ~a~%" index file-name))
          do (incf (total-sequence-count *original-analysis*))
          when (and num (>= index num))
            do (return-from read-loop nil)
          if (accept-result res)
            do (let ((key (mapcar (lambda (digit) (del-name (row-code digit))) res)))
                 (incf (good-sequence-count *original-analysis*))
                 (incf (gethash key (sequences *original-analysis*) 0))
                 (when pass-file (sa:write-record pass-file title seq)))
          else
            do (when fail-file (sa:write-record fail-file title seq)))
    (format t "Read ~a total sequences and ~a good sequences - there are ~a unique sequences~%"
            num (good-sequence-count *original-analysis*)
            (hash-table-count (sequences *original-analysis*)))))

(defun analyze (seq-file-path &key (num 10) verbose pass-file-name fail-file-name)
  (unless (probe-file seq-file-path)
    (error "Could not find file ~a" seq-file-path))
  (let ((seq-file (sa:make-seq-file-in seq-file-path)))
    (let ((pass-file (when pass-file-name
                       (sa:make-seq-file-out pass-file-name)))
          (fail-file (when fail-file-name
                       (sa:make-seq-file-out fail-file-name))))
      (analyze-seq-file seq-file-path seq-file :num num :verbose verbose
                                               :pass-file pass-file
                                               :fail-file fail-file))))

(defun serial-analyze (sequence-files &key (num 10) pass-file-name fail-file-name)
  (clear-analysis sequence-files)
  (loop for file in sequence-files
        do (unless (probe-file file)
             (error "Could not find file: ~a~%" file)))
  (loop for file in sequence-files
        do (analyze file :num num :pass-file-name pass-file-name :fail-file-name fail-file-name))
  (format t "Finished.~%"))
  

(defun parallel-analyze (sequence-files &key (num 10))
  (clear-analysis)
  (loop for file in sequence-files
        do (unless (probe-file file)
             (error "Could not find file: ~a~%" file)))
  (let ((threads (loop for tinum from 0
                       for file in sequence-files
                       for thread = (mp:process-run-function
                                      (format nil "process~d" tinum)
                                      (lambda ()
                                        (analyze file :num num)))
                       collect thread)))
    (loop for thread in threads
          do (mp:process-join thread))))

(defun save-csv (vals file)
  (with-open-file (fout file :direction :output)
    (loop for row in vals
          do (let ((*print-base* 16))
               (format fout "~{ ~s~^,~}," (car row)))
          do (format fout "~{ ~s~^,~}~%" (cdr row)))))

(defun sequence-survey (seq &optional (analysis *original-analysis*))
  (let ((ht (make-hash-table :test 'equal)))
    (maphash (lambda (k v)
               (let ((seq-key (subseq k 1 4)))
                 (push (cons k v) (gethash seq-key ht))))
             analysis)
    (gethash seq ht)))

(defun show-hash-table (result &optional (max-show 16))
  (let ((count max-show))
    (block showloop
      (maphash (lambda (k v)
                 (when (zerop count)
                   (return-from showloop nil))
                 (format t "~s ~a~%" k v)
                 (decf count))
               result)))
  (format t "There are ~a entries~%" (hash-table-count result)))


(defun search-library-code (library-code sequences &key (max-show 20))
  (let ((result (make-hash-table :test 'equal)))
    (maphash (lambda (k v)
               (let ((lib-key (elt k 4)))
                 (when (string= lib-key library-code)
                   (setf (gethash k result) v))))
             sequences)
    (show-hash-table result max-show)
    result))

(defun library-codes (sequences)
  "Gather up all of the unique library codes and count how many times they appear"
  (let ((library-codes (make-hash-table :test 'equal)))
    (maphash (lambda (key value)
               (let ((library-code (subseq key 4 5)))
                 (incf (gethash library-code library-codes 0) value)))
             sequences)
    (maphash (lambda (key value)
               (format t "~s ~a~%" key value))
             library-codes)
    nil))

(defun bead-codes ()
  "Gather the sequences and count how many times they appear with different bead codes"
  (let ((bead-codes (make-hash-table :test 'equal)))
    (maphash (lambda (key value)
               (let ((library-code (subseq key 0 2)))
                 (incf (gethash library-code bead-codes 0) value)))
             (sequences *original-analysis*))
    bead-codes))
    
(defun save-results (file &key force)
  "Save the results to a file. Pass :force t if you want to force overwriting an existing file"
  (when (and (not force) (probe-file file))
    (error "The file ~a already exists~%" file))
  (with-open-file (fout file :direction :output)
    (let ((*print-readably* t)
          (*print-pretty* nil))
      (princ (file-names *original-analysis*) fout)
      (terpri fout)
      (princ (total-sequence-count *original-analysis*) fout)
      (terpri fout)
      (princ (good-sequence-count *original-analysis*) fout)
      (terpri fout)
      (princ (hash-table-count (sequences *original-analysis*)) fout)
      (terpri fout)
      (maphash (lambda (k v)
                 (format fout "~s ~s~%" k v))
               (sequences *original-analysis*))))
  (format t "Saved to ~a~%" file))

(defun load-results (file)
  (with-open-file (fin file :direction :input)
    (let ((analysis (make-instance 'analysis)))
      (setf (file-names analysis) (read fin))
      (setf (total-sequence-count analysis) (read fin))
      (setf (good-sequence-count analysis) (read fin))
      (let* ((num-pairs (read fin)))
        (loop for index below num-pairs
              for key = (read fin)
              for val = (read fin)
              do (setf (gethash key (sequences analysis)) val)))
      analysis)))

(defun join-codons (input)
  (let ((analysis (make-hash-table :test 'equal)))
    (maphash (lambda (key val)
               (let* ((sequence-key (subseq key 0 10))
                      (compact-sequence-key (loop for cur = sequence-key then (cddr cur)
                                                  for msd = (car cur)
                                                  for lsd = (cadr cur)
                                                  when (null cur)
                                                    return digits
                                                  when cur
                                                    collect (join-del-codes msd lsd) into digits
                                                  )))
                 (setf (gethash compact-sequence-key analysis) val)))
             input)
    analysis))

(defun sequence-counts (input)
  (let ((counts nil))
    (maphash (lambda (k v)
               (declare (ignore k))
               (push (log v 10.0) counts))
             input)
    (let ((sorted-counts (sort counts #'<)))
      (loop for index from 0
            for count in sorted-counts
            collect (cons index count)))))

(defun filter-sequences (input &optional (threshold *minimum-counts*))
  (let ((result-sequences (make-hash-table :test 'equal)))
    (maphash (lambda (key value)
               (when (> value threshold)
                 (setf (gethash key result-sequences) value)))
             input)
    result-sequences))

(defun bead-specific-counts (input-analysis &optional (key-positions '(1 2 3)))
  (let ((analysis (make-hash-table :test 'equal)))
    (maphash (lambda (key val)
               (declare (ignore val))
               (let ((sequence-key (loop for pos in key-positions
                                         collect (elt key pos)))) 
                 (incf (gethash sequence-key analysis 0))))
             input-analysis)
    analysis))

(defun merge-hits (in1 in2)
  "Merge hits from two sorts using their sequences"
  (let ((results nil))
    (maphash (lambda (k1 v1)
               (let ((v2 (gethash k1 in2 0)))
                 (push (list k1 v1 v2) results)))
             in1)
    (sort results (lambda (x y)
                    (if (> (second x) (second y))
                        t
                        (if (= (second x) (second y))
                            (> (third x) (third y))))))))

(defvar *hits* nil)
(defun sort-hits (&optional analysis)
  (let ((unsorted))
    (maphash (lambda (key value)
               (push (cons key value) unsorted))
             analysis)
    (setf *hits* (sort unsorted #'> :key #'cdr))
    *hits*))

(defun multiple-redundant (hits &key (min 2))
  "Return sequences that are multiply redundant.
Arguments:
hits - A list of merged hits.
min  - The minimum number of redundancies."
  (loop for row in hits
        for red1 = (second row)
        for red2 = (third row)
        when (and (>= red1 min)
                  (>= red2 min))
          collect row))
