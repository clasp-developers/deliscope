(in-package :deliscope)

(defparameter *simple-score* (sa:make-simple-score 0 -1 -2))
(defparameter *align-config* (sa:make-align-config :t-nil-nil-t))
(defparameter *del-code-format* "~d~a~2,'0d")

(defparameter *minimum-phred-quality* 20)
(defparameter *low-quality-max* 2)
(defparameter *minimum-counts* 100)

(defclass align-view ()
  ((align :initarg :align :accessor align)
   (row0& :initarg :row0& :reader row0&)
   (row1& :initarg :row1& :reader row1&)))
   
(defclass parser ()
  ((name :initarg :name :accessor name)
   (files :initarg :files :accessor files)
   (num-sequences-per-file :initarg :num-sequences-per-file :accessor num-sequences-per-file)
   (pass-file-name :initarg :pass-file-name :accessor pass-file-name)
   (fail-file-name :initarg :fail-file-name :accessor fail-file-name)
   (output-file-name :initarg :output-file-name :accessor output-file-name)
   (overwrite :initarg :overwrite :accessor overwrite)
   ))

(defclass filtering ()
  ((name :initarg :name :accessor name)
   (file-names :initarg :file-names :accessor file-names)
   (sequences :initform (make-hash-table :test 'equal) :initarg :sequences :accessor sequences)
   (start-time :initarg :start-time :accessor start-time)
   (current-time :initarg :current-time :accessor current-time)
   (total-sequence-count :initform 0 :initarg :total-sequence-count :accessor total-sequence-count)
   (good-sequence-count :initform 0 :initarg :good-sequence-count :accessor good-sequence-count)))

(defmethod summarize (filtering)
  (format t "Files:~%~s~%" (file-names filtering))
  (format t "Total sequences: ~a~%" (total-sequence-count filtering))
  (format t "Good sequences:  ~a~%" (good-sequence-count filtering))
  (format t "Sample sequences:~%")
  (show-hash-table (sequences filtering))
  (values))

(defclass job-tracker ()
  ((sequence-count :initform 0 :accessor sequence-count)
   (good-sequence-count :initform 0 :accessor good-sequence-count)
   (seq-file-position :initform 0 :accessor seq-file-position)
   (seq-file-length :initform 0 :accessor seq-file-length)
   (result :accessor result)))

(defun sum-trackers (sum-tracker trackers)
  (loop for index below (length trackers)
        for tracker = (elt trackers index)
        sum (sequence-count tracker) into sequence-count-sum
        sum (good-sequence-count tracker) into good-sequence-count-sum
        sum (seq-file-position tracker) into seq-file-position-sum
        finally (setf (sequence-count sum-tracker) sequence-count-sum
                      (good-sequence-count sum-tracker) good-sequence-count-sum
                      (seq-file-position sum-tracker) seq-file-position-sum)))


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
  ((score :initform nil :initarg :score :accessor score)
   (low-quality-count :initform nil :initarg :low-quality-count :accessor low-quality-count)
   (quality-string :initform nil :initarg :quality-string :accessor quality-string)
   (row-code :initform nil :initarg :row-code :accessor row-code)))

(defmethod print-object ((object match) stream)
  (print-unreadable-object (object stream :type t)
    (format stream "~a ~a ~a ~s"
            (score object)
            (row-code object)
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
    (princ (with-output-to-string (sout)
             (let ((*print-pretty* t))
               (format sout "---- 1x sequences hamming distances~%~s~%" (hamming-matrix r1xseqs))
               (format sout "---- 2x sequences hamming distances~%~a~%" (hamming-matrix r2xseqs)))))
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
      (format t "---- Relative codon offsets:~% ~{~a ~}~%" (mapcar (lambda (cd) (format nil "~a...~a" (start cd) (end cd))) column-codes))
      (setf *all-codes* (make-instance 'all-codes
                                       :forward-primer (sa:make-string :dna5q-string forward-primer)
                                       :column-codes column-codes
                                       :end (end (car (last column-codes))))))))


(defparameter *library-code* nil)

(defun set-filter (&key library-code (minimum-counts 100))
  "Set the library code to use as a filter for sequences"
  (setf *library-code* library-code)
  (setf *minimum-counts* minimum-counts))

(defun analyze-column-using-align (align-view sequence column-codes start adjust low-quality-count quality-string &key verbose)
  (let* ((abs-start (+ start (start column-codes) adjust))
         (best-score -999999)
         (best-row-code nil))
    (loop for row-code in (codes column-codes)
          for name = (del-name row-code)
          for code = (code row-code)
          for lower-diag = (- abs-start 1)
          for upper-diag = (+ lower-diag 1) ;;+ start end-offset )
          for score = (progn
                        (sa:assign-source (row0& align-view) sequence)
                        (sa:assign-source (row1& align-view) code)
                        (sa:global-alignment (align align-view) *simple-score* *align-config*
                                                          lower-diag
                                                          upper-diag ))
          do (when verbose
               (format t "Score: ~a~%" score)
               (format t "~a~%" (sa:to-string (align align-view))))
          do (when (> score best-score)
               (setf best-score score
                     best-row-code row-code)))
    (make-instance 'match
                   :score best-score
                   :low-quality-count low-quality-count
                   :quality-string quality-string
                   :row-code best-row-code)))

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

(defun analyze-column (align-view sequence column-codes start adjust &key verbose)
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
      (analyze-column-using-align align-view sequence column-codes start adjust low-quality-count qual-string :verbose verbose))))

(defun analyze-forward-primer (align-view sequence &key verbose)
  (sa:assign-source (row0& align-view) sequence)
  (sa:assign-source (row1& align-view) (forward-primer *all-codes*))
  (let ((fwd-score (sa:global-alignment (align align-view) *simple-score* *align-config*))
        (start (1+ (sa:to-view-position (row1& align-view) (sa:length (forward-primer *all-codes*))))))
    (when verbose (format t "Alignment: ~%~a~%" (sa:to-string (align align-view))))
    (values start fwd-score)))

(defun analyze-sequence (align-view sequence &key verbose)
  (multiple-value-bind (start fwd-score)
      (analyze-forward-primer align-view sequence :verbose verbose)
    (declare (ignore fwd-score))
    (when verbose
      (format t "len fwd: ~a~%" (sa:length (forward-primer *all-codes*)))
      (format t "3' end: ~a~%" start)
      (format t "~a~%" (sa:to-string (align align-view))))
    (when (> (sa:length sequence) (+ start (end *all-codes*)))
      (let ((result (loop for column-codes in (column-codes *all-codes*)
                          for start-offset = (start column-codes)
                          for end-offset = (end column-codes)
                          for best-digit = (analyze-column align-view sequence column-codes start 0 :verbose verbose)
                          do (when verbose (format t "offset: ~a best-digits: ~a~%" start-offset best-digit))
                          collect best-digit)))
        (when verbose
          (format t "result: ~a~%" result))
        result))))


(defvar *filtering* (make-hash-table :test 'equal))
(defvar *bead-specific (make-hash-table :test 'equal))


(defun accept-result (result)
  (and result
       (every (lambda (x) (>= (score x) -2)) result)
       (every (lambda (x) (< (low-quality-count x) *low-quality-max*)) result)))


(defun evaluate-progress (start-time sum-tracker total-count)
  (let ((count (sequence-count sum-tracker)))
    (if (= count 0)
        (let ((msg (format nil "Read ~a of ~a~%Estimating remaining time...~%" count total-count)))
          (values 0 msg))
        (let* ((good-sequence-count (good-sequence-count sum-tracker))
               (current-time (get-internal-real-time))
               (elapsed-seconds (float (/ (- current-time start-time) internal-time-units-per-second)))
               (remaining-seconds (if (= count 0) 9999999.0 (- (* (/ total-count count) elapsed-seconds) elapsed-seconds)))
               (remaining-minutes (/ remaining-seconds 60.0))
               (remaining-hours (/ remaining-minutes 60.0)))
          (multiple-value-bind (hours minutes)
              (if (< remaining-hours 0.0)
                  (values 0 0)
                  (floor remaining-hours))
            (labels ((time-description (hours minutes)
                       (if (zerop hours)
                           (format nil "~4,1f minute~:p" minutes)
                           (format nil "~a hour~:p ~4,1f minute~:p" hours minutes))))
              (let ((msg (format nil "Read ~a of ~a~%Good sequences: ~a~%Estimated remaining time: ~a~%"
                                 count total-count
                                 good-sequence-count
                                 (time-description hours (* 60.0 minutes))))
                    (progress (floor (float (* 100.0 (/ count total-count))))))
                (values progress msg))))))))

(defparameter *align-lock* (bt:make-lock 'align-lock))
(defun thread-safe-make-align-view (type)
  (unwind-protect
       (progn
         (bt:acquire-lock *align-lock*)
         (let* ((align (sa:make-align type))
                (rows& (sa:rows& align)))
           (sa:resize rows& 2)
           (let ((row0& (sa:row& align 0))
                 (row1& (sa:row& align 1)))
             (make-instance 'align-view
                            :align align
                            :row0& row0&
                            :row1& row1&))))
    (bt:release-lock *align-lock*)))

(defparameter *string-lock* (bt:make-lock 'string-lock))
(defun thread-safe-make-string (type)
  (unwind-protect
       (progn
         (bt:acquire-lock *string-lock*)
         (sa:make-string type))
    (bt:release-lock *string-lock*)))

(defun analyze-one-file (file-index file-name local-file-count &key progress-callback pass-fail-files verbose)
  (let ((sequences (make-hash-table :test 'equal))
        (seq-file (sa:make-seq-file-in file-name))
        (pass-file nil)
        (fail-file nil)
        (count 0)
        (good-sequence-count 0))
    (when pass-fail-files
      (let* ((pass-file-name (make-pathname
                              :name (format nil "~a-pass" (pathname-name (pathname file-name)))
                              :type "fastq" :defaults (pathname file-name)))
             (fail-file-name (make-pathname
                              :name (format nil "~a-fail" (pathname-name (pathname file-name)))
                              :type "fastq" :defaults (pathname file-name))))
        (setf pass-file (sa:make-seq-file-out (namestring pass-file-name))
              fail-file (sa:make-seq-file-out (namestring fail-file-name)))))
    (unwind-protect
         (let ((title (thread-safe-make-string :char-string))
               (seq (thread-safe-make-string :dna5q-string))
               (align-view (thread-safe-make-align-view :dna5q-string))
               (next-progress-update (if (> (+ count 10000) local-file-count)
                                         local-file-count
                                         (+ count 10000))))
           (loop named read-loop
                 for index from 0
                 for res = (progn
                             (when (sa:at-end seq-file)
                               (return-from read-loop count))
                             (sa:read-record title seq seq-file)
                             (analyze-sequence align-view seq :verbose verbose))
                 when (>= count next-progress-update)
                   do (progn
                        (incf next-progress-update 10000)
                        (when (> next-progress-update local-file-count)
                          (setf next-progress-update local-file-count))
                        (when progress-callback
                          (funcall progress-callback file-index count good-sequence-count 9999)))
                 do (incf count)
                 unless (< index local-file-count)
                   do (return-from read-loop count)
                 if (accept-result res)
                   do (let ((key (mapcar (lambda (digit) (del-name (row-code digit))) res)))
                        (incf good-sequence-count)
                        (incf (gethash key sequences 0))
                        (when pass-file (sa:write-record pass-file title seq)))
                 else
                   do (when fail-file (sa:write-record fail-file title seq))
                 finally (return-from read-loop count)))
      (progn
        (sa:close seq-file)
        (when pass-file (sa:close pass-file))
        (when fail-file (sa:close fail-file))))
    sequences))


(defun make-job-tracker-vector (number-of-jobs)
  (let ((trackers (make-array number-of-jobs)))
    (loop for index below number-of-jobs
          do (setf (elt trackers index) (make-instance 'job-tracker)))
    trackers))

(defun serial-analyze (parser &key progress-callback)
  "Read sequences and generate an filtering"
  (let* ((start-time (get-internal-real-time))
         (filtering (make-instance 'filtering
                                   :name (name parser)
                                   :start-time start-time
                                   :file-names (files parser))))
    (loop for file in (files parser)
          do (unless (probe-file file)
               (error "Could not find file: ~a~%" file)))
    (let ((total-count (apply '+ (num-sequences-per-file parser)))
          (trackers (make-job-tracker-vector (length (files parser))))
          (sum-tracker (make-instance 'job-tracker)))
      (loop for file in (files parser)
            for file-index from 0
            for local-file-count in (num-sequences-per-file parser)
            for one-file-ht = (analyze-one-file
                               file-index file local-file-count
                               :progress-callback
                               (lambda (file-index count good-sequence-count seq-file-position)
                                 (let ((tracker (elt trackers file-index)))
                                   (setf (sequence-count tracker) count
                                         (good-sequence-count tracker) good-sequence-count
                                         (seq-file-position tracker) seq-file-position)
                                   (when progress-callback
                                     (sum-trackers sum-tracker trackers)
                                     (multiple-value-bind (progress msg)
                                         (evaluate-progress start-time sum-tracker total-count)
                                       (funcall progress-callback progress msg)))))
                               :pass-fail-files t)
            do (maphash (lambda (seq counts)
                          (incf (gethash seq (sequences filtering) 0) counts))
                        one-file-ht)))
    (format t "About to save results~%")
    (save-filtering filtering (output-file-name parser) :overwrite (overwrite parser))
    (when progress-callback (funcall progress-callback 100 nil :done t))
    (format t "Finished filtering.~%")))


(defparameter *progress-callback* nil)
(defparameter *file-index* nil)
(defparameter *file* nil)
(defparameter *num-sequences* nil)
(defparameter *trackers* nil)

(defun parallel-analyze (parser &key progress-callback)
  "Read sequences and generate an filtering"
  (let* ((start-time (get-internal-real-time))
         (filtering (make-instance 'filtering
                                   :name (name parser)
                                   :start-time start-time
                                   :file-names (files parser))))
    (loop for file in (files parser)
          do (unless (probe-file file)
               (error "Could not find file: ~a~%" file)))
    (let* ((trackers (make-job-tracker-vector (length (files parser))))
           (sum-total (apply #'+ (num-sequences-per-file parser)))
           (workers (loop for file-index below (length (files parser))
                          for file in (files parser)
                          for num-sequences-in-file in (num-sequences-per-file parser)
                          collect (bt:make-thread
                                   (lambda ()
                                     (let ((ht (analyze-one-file
                                                *file-index* *file* *num-sequences*
                                                :progress-callback (lambda (file-index count good-sequence-count seq-file-position)
                                                                     (let ((tracker (elt *trackers* file-index)))
                                                                       (setf (sequence-count tracker) count
                                                                             (good-sequence-count tracker) good-sequence-count
                                                                             (seq-file-position tracker) seq-file-position))))))
                                       (setf (result (elt *trackers* *file-index*)) ht)))
                                   :initial-bindings (list (cons '*file-index* file-index)
                                                           (cons '*file* file)
                                                           (cons '*num-sequences* num-sequences-in-file)
                                                           (cons '*trackers* trackers))))))
      (let ((sum-tracker (make-instance 'job-tracker)))
        (loop named monitor
              do (sleep 1)
                 (core:check-pending-interrupts)
                 (sum-trackers sum-tracker trackers)
                 (when progress-callback
                   (multiple-value-bind (progress msg)
                       (evaluate-progress start-time sum-tracker sum-total)
                     (funcall progress-callback progress msg)))
                 (when (every (lambda (x) (not (bt:thread-alive-p x))) workers)
                   (return-from monitor nil))))
      (loop for tracker across trackers
            for result = (result tracker)
            do (maphash (lambda (seq counts)
                          (incf (gethash seq (sequences filtering) 0) counts))
                        result))
      (format t "About to save results~%")
      (save-filtering filtering (output-file-name parser) :overwrite (overwrite parser))
      (when progress-callback (funcall progress-callback 100 nil :done t))
      (format t "Finished parsing sequences.~%"))))

(defun sequence-survey (seq filtering)
  (let ((ht (make-hash-table :test 'equal)))
    (maphash (lambda (k v)
               (let ((seq-key (subseq k 1 4)))
                 (push (cons k v) (gethash seq-key ht))))
             filtering)
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

(defun search-library-codes (sequences)
  "Gather up all of the unique library codes and count how many times they appear"
  (let ((library-codes (make-hash-table :test 'equal)))
    (maphash (lambda (key value)
               (let ((library-code (subseq key 4 5)))
                 (incf (gethash (first library-code) library-codes 0) value)))
             sequences)
    library-codes))

(defun bead-codes (filtering)
  "Gather the sequences and count how many times they appear with different bead codes"
  (let ((bead-codes (make-hash-table :test 'equal)))
    (maphash (lambda (key value)
               (let ((library-code (subseq key 0 2)))
                 (incf (gethash library-code bead-codes 0) value)))
             (sequences filtering))
    bead-codes))
    
(defun save-filtering (filtering file &key overwrite)
  "Save the results to a file. Pass :force t if you want to force overwriting an existing file"
  (when (and (not overwrite) (probe-file file))
    (error "The file ~a already exists~%" file))
  (with-open-file (fout file :direction :output :if-exists (if overwrite :supersede))
    (let ((*print-readably* t)
          (*print-pretty* nil))
      (prin1 (name filtering) fout)
      (terpri fout)
      (prin1 (file-names filtering) fout)
      (terpri fout)
      (prin1 (total-sequence-count filtering) fout)
      (terpri fout)
      (prin1 (good-sequence-count filtering) fout)
      (terpri fout)
      (prin1 (hash-table-count (sequences filtering)) fout)
      (terpri fout)
      (maphash (lambda (k v)
                 (format fout "~s ~s~%" k v))
               (sequences filtering))))
  (format t "Saved to ~a~%" file))

(defun load-filtering (file)
  (with-open-file (fin (merge-pathnames file) :direction :input)
    (let ((filtering (make-instance 'filtering)))
      (setf (name filtering) (read fin))
      (setf (file-names filtering) (read fin))
      (setf (total-sequence-count filtering) (read fin))
      (setf (good-sequence-count filtering) (read fin))
      (let* ((num-pairs (read fin)))
        (loop for index below num-pairs
              for key = (read fin)
              for val = (read fin)
              do (setf (gethash key (sequences filtering)) val)))
      filtering)))

(defun join-codons (input)
  (let ((filtering (make-hash-table :test 'equal)))
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
                 (setf (gethash compact-sequence-key filtering) val)))
             input)
    filtering))

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

(defun bead-specific-counts (input-filtering &key (key-positions '(1 2 3)) library-code)
  (let ((filtering (make-hash-table :test 'equal)))
    (maphash (lambda (key val)
               (declare (ignore val))
               (let ((library-key (elt key 4))
                     (sequence-key (loop for pos in key-positions
                                         collect (elt key pos)))) 
                 (if library-code
                     (when (string= library-code library-key)
                       (incf (gethash sequence-key filtering 0)))
                     (incf (gethash sequence-key filtering 0)))))
             input-filtering)
    filtering))

(defclass analysis (filtering)
  ((joined-sequences :initarg :jointed-sequences :accessor joined-sequences)
   (minimum-counts :initform *minimum-counts* :accessor minimum-counts)
   (filtered-sequences :accessor filtered-sequences)
   (library-code :initform *library-code* :accessor library-code)
   (keys :initarg :keys :accessor keys)
   (redundancy :accessor redundancy)
   ))

(defun process-results (file &key (keys '(1 2 3)))
  (let* ((analysis (load-filtering file)))
    (change-class analysis 'analysis :keys keys)
    (setf (joined-sequences analysis) (join-codons (sequences analysis)))
    (setf (filtered-sequences analysis) (filter-sequences (joined-sequences analysis) (minimum-counts analysis)))
    (setf (redundancy analysis) (bead-specific-counts (filtered-sequences analysis)
                                                      :key-positions keys
                                                      :library-code (library-code analysis)))
    analysis))


(defun filtered-library-codes (analysis)
  (let ((lib-codes (search-library-codes (filtered-sequences analysis))))
    (maphash (lambda (k v)
               (format t "~s ~a~%" k v))
             lib-codes))
  (values))



(defun compare-list (minimum-redundancy in1 &optional in2)
  "Compare hits from two sorts using their sequences"
  (let ((results nil))
    (maphash (lambda (k1 v1)
               (when (>= v1 minimum-redundancy)
                 (if in2
                     (let ((v2 (gethash k1 (redundancy in2) 0)))
                       (when (> v2 minimum-redundancy)
                         (push (list k1 v1 v2) results)))
                     (push (list k1 v1) results))))
             (redundancy in1))
    (if in2
        (sort results (lambda (x y)
                        (if (> (second x) (second y))
                            t
                            (if (= (second x) (second y))
                                (> (third x) (third y))))))
        (sort results (lambda (x y)
                        (> (second x) (second y)))))))

(defun compare (in1 in2 &key (minimum-redundancy 1))
  (compare-list minimum-redundancy in1 in2))

(defun redundancies (in1 &key (minimum-redundancy 1))
  (compare-list minimum-redundancy in1))

(defvar *hits* nil)
(defun sort-hits (&optional filtering)
  (let ((unsorted))
    (maphash (lambda (key value)
               (push (cons key value) unsorted))
             filtering)
    (setf *hits* (sort unsorted #'> :key #'cdr))
    *hits*))

(defun save-csv (file analysis1 &optional analysis2)
  (with-open-file (fout file :direction :output)
    (format fout "~{ ~a~^,~}~%" (list* "Pos1" "Pos2" "Pos3"
                                       (format nil "~a" (name analysis1))
                                       (if analysis2
                                           (list (format nil "~a" (name analysis2)))
                                           nil)))
    (let ((rows (compare-list 0 analysis1 analysis2)))
      (loop for row in rows
            do (format fout "~{ ~a~^,~}," (car row))
            do (format fout "~{ ~a~^,~}~%" (cdr row)))))
  (format t "Wrote to ~a~%" file)
  (values))


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


#|

Copyright (c) 2021, Christian E. Schafmeister

Deliscope is free software; you can redistribute it and/or
modify it under the terms of the GNU Library General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

See file 'LICENSE' for full details.

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

|#
