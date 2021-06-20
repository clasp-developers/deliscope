(in-package :deliscope)

(defparameter *simple-score* (sa:make-simple-score 0 -1 -2))
(defparameter *align-config* (sa:make-align-config :t-nil-nil-t))
(defparameter *del-code-format* "~d~a~2,'0d")

(defparameter *minimum-phred-quality* 20)
(defparameter *low-quality-max* 2)
(defparameter *minimum-counts* 100)

(defclass outputs ()
  ((pass-file :initarg :pass-file :accessor pass-file)
   (pass-file-lock :initarg :pass-file-lock :accessor pass-file-lock)
   (fail-file :initarg :fail-file :accessor fail-file)
   (fail-file-lock :initarg :fail-file-lock :accessor fail-file-lock)))

(defclass sequence-align ()
  ((metadata :initarg :metadata :accessor metadata)
   (seq :initarg :seq :accessor seq)
   (align :initarg :align :accessor align)
   (row0& :initarg :row0& :reader row0&)
   (row1& :initarg :row1& :reader row1&)
   (pass-fail :accessor pass-fail)
   (results-lock :accessor results-lock)
   (results :accessor results)))

(defparameter *seqan-make-lock* (bt:make-lock 'seqan-make-lock))
(defun make-sequence-align ()
  "seqan doesn't appear to have thread safe allocators - so we use a lock
when allocating align and string objects"
  (unwind-protect
       (progn
         (bt:acquire-lock *seqan-make-lock*)
         (let ((align (seqan:make-align :dna5q-string)))
           (seqan:resize (sa:rows& align) 2)
           (make-instance 'sequence-align
                          :metadata (seqan:make-string :char-string)
                          :seq (seqan:make-string :dna5q-string)
                          :align align
                          :row0& (seqan:row& align 0)
                          :row1& (seqan:row& align 1))))
    (bt:release-lock *seqan-make-lock*)))
                 
(defclass parser ()
  ((name :initarg :name :accessor name)
   (files :initarg :files :accessor files)
   (num-sequences-per-file :initarg :num-sequences-per-file :accessor num-sequences-per-file)
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

(defclass worker ()
  ((resource-stack :initarg :resource-stack :reader resource-stack)
   (tracker :initarg :tracker :reader tracker)
   (worker-index :initarg :worker-index :reader worker-index)
   (thread :initform nil :accessor thread)))

(defclass job-tracker ()
  ((sequence-count :initform 0 :accessor sequence-count)
   (good-sequence-count :initform 0 :accessor good-sequence-count)))

(defun sum-trackers (sum-tracker trackers)
  (loop for index below (length trackers)
        for tracker = (elt trackers index)
        sum (sequence-count tracker) into sequence-count-sum
        sum (good-sequence-count tracker) into good-sequence-count-sum
        finally (setf (sequence-count sum-tracker) sequence-count-sum
                      (good-sequence-count sum-tracker) good-sequence-count-sum)))

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

(defmacro vformat (verbose stream fmt &rest args)
  `(when ,verbose
     (bt:with-lock-held (,verbose)
       (format ,stream ,fmt ,@args))))

(defun analyze-column-using-align (sequence-align sequence column-codes start adjust low-quality-count quality-string &key verbose)
  (let* ((abs-start (+ start (start column-codes) adjust))
         (best-score -999999)
         (best-row-code nil))
    (loop for row-code in (codes column-codes)
          for name = (del-name row-code)
          for code = (code row-code)
          for lower-diag = (- abs-start 1)
          for upper-diag = (+ lower-diag 1) ;;+ start end-offset )
          for score = (progn
                        (sa:assign-source (row0& sequence-align) sequence)
                        (sa:assign-source (row1& sequence-align) code)
                        (sa:global-alignment (align sequence-align) *simple-score* *align-config*
                                                          lower-diag
                                                          upper-diag ))
          do (vformat verbose t "Score: ~a~%" score)
             (vformat verbose t "~a~%" (sa:to-string (align sequence-align)))
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
    (vformat verbose t "abs-start: ~d  abs-end ~d Subsequence: ~a~%" abs-start abs-end (subseq (sa:to-string sequence) (1- abs-start) (1- abs-end)))
    (loop for row-code in (codes column-codes)
          for code-seq = (sequence row-code)
          when (string= seq-string code-seq :start1 (1- abs-start) :end1 (1- abs-end))
            do (return-from analyze-column-fast (values t row-code)))
    (values nil)))

(defun analyze-column (sequence-align sequence column-codes start adjust &key verbose)
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
      (analyze-column-using-align sequence-align sequence column-codes start adjust low-quality-count qual-string :verbose verbose))))

(defun analyze-forward-primer (sequence-align sequence &key verbose)
  (sa:assign-source (row0& sequence-align) sequence)
  (sa:assign-source (row1& sequence-align) (forward-primer *all-codes*))
  (let ((fwd-score (sa:global-alignment (align sequence-align) *simple-score* *align-config*))
        (start (1+ (sa:to-view-position (row1& sequence-align) (sa:length (forward-primer *all-codes*))))))
    (vformat verbose t "Alignment: ~%~a~%" (sa:to-string (align sequence-align)))
    (values start fwd-score)))

(defun analyze-sequence (sequence-align &key verbose worker)
  (vformat verbose t "worker: ~a analyze-sequence for ~a~%" worker sequence-align)
  (let ((sequence (seq sequence-align)))
    (multiple-value-bind (start fwd-score)
        (analyze-forward-primer sequence-align sequence :verbose verbose)
      (declare (ignore fwd-score))
      (vformat verbose t "worker: ~a len fwd: ~a~%" worker (sa:length (forward-primer *all-codes*)))
      (vformat verbose t "worker: ~a 3' end: ~a~%" worker start)
      (vformat verbose t "worker: ~a ~a~%" worker (sa:to-string (align sequence-align)))
      (when (> (sa:length sequence) (+ start (end *all-codes*)))
        (let ((result (loop for column-codes in (column-codes *all-codes*)
                            for start-offset = (start column-codes)
                            for end-offset = (end column-codes)
                            for best-digit = (analyze-column sequence-align sequence column-codes start 0 :verbose verbose)
                            do (vformat verbose t "worker: ~a offset: ~a best-digits: ~a~%" worker start-offset best-digit)
                            collect best-digit)))
          (vformat verbose t "worker: ~a result: ~a~%" worker result)
        result)))))


(defvar *filtering* (make-hash-table :test 'equal))
(defvar *bead-specific (make-hash-table :test 'equal))


(defun accept-result (result)
  (and result
       (every (lambda (x) (>= (score x) -2)) result)
       (every (lambda (x) (< (low-quality-count x) *low-quality-max*)) result)))


(defun evaluate-progress (start-time count good-sequence-count total-count file-sequence-count file)
  (if (= count 0)
      (let ((msg (format nil "Read: ~a from ~a~%Parsed ~a of ~a~%Estimating remaining time...~%" file-sequence-count file count total-count)))
        (values 0 msg))
      (let* ((current-time (get-internal-real-time))
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
            (let ((msg (format nil "Read ~a sequences from ~a~%Parsed ~a of ~a total sequences~%Good sequences: ~a~%Estimated remaining time: ~a~%"
                               file-sequence-count file
                               count total-count
                               good-sequence-count
                               (time-description hours (* 60.0 minutes))))
                  (progress (floor (float (* 100.0 (/ count total-count))))))
              (values progress msg)))))))


(defun wait-for-sequence (worker queue &key verbose)
  (let ((tracker (tracker worker)))
    (symbol-macrolet ((sequence-count (sequence-count tracker))
                      (good-sequence-count (good-sequence-count tracker)))
      (loop named work-loop
            for sequence-align = (let ((sal (progn
                                              (vformat verbose t "worker: ~a About to pop queue~%" (worker-index worker))
                                              (lparallel.queue:pop-queue queue))))
                                   (when (eq sal :quit)
                                     (vformat verbose t "worker ~a exiting~%" (worker-index worker))
                                     (return-from work-loop nil))
                                   (unless sal
                                     (vformat verbose t "worker: ~a got nil as sal~%" (worker-index worker)))
                                   sal)
            for metadata = (metadata sequence-align)
            for seq = (seq sequence-align)
            for pass-fail = (pass-fail sequence-align)
            for result = (analyze-sequence sequence-align :verbose nil :worker (worker-index worker))
            do (mp:atomic-incf (slot-value tracker 'sequence-count))
            if (and result (accept-result result))
              do (let ((key (mapcar (lambda (digit) (del-name (row-code digit))) result)))
                   (mp:atomic-incf (slot-value tracker 'good-sequence-count))
                   (bt:with-lock-held ((results-lock sequence-align))
                     (vformat verbose t "worker: ~a recording hit: ~a~%" (worker-index worker) key)
                     (incf (gethash key (results sequence-align) 0))
                     (vformat verbose t "worker: ~a results = ~s~%" (worker-index worker) (results sequence-align)))
                   (when pass-fail
                     (mp:with-lock ((pass-file-lock pass-fail))
                       (sa:write-record (pass-file pass-fail) metadata seq))))
            else
              do (when pass-fail
                   (mp:with-lock ((fail-file-lock pass-fail))
                     (sa:write-record (fail-file pass-fail) metadata seq)))
            do (vformat verbose t "worker: ~a atomic-push of resource back on stack~%" (worker-index worker))
            do (mp:atomic-push sequence-align (car (resource-stack worker))) ; recycle the sequence-align 
            ))))


(defparameter *worker* nil)
(defparameter *queue* nil)
(defparameter *verbose* nil)
(defun analyze-parsers-using-workers (number-of-workers parsers
                                      &key max-sequences-per-file progress-callback pass-fail-files verbose testing)
  (vformat verbose t "Entered analyze-files-using-workers ~a~%" parsers)
  ;; Create resources and workers
  (format t "max-sequences-per-file: ~a~%" max-sequences-per-file)
  (let* ((sequence-count 0)
         (start-time (get-internal-real-time))
         (estimated-sequences (if max-sequences-per-file
                                  (loop for parser in parsers
                                        sum (loop for num in (num-sequences-per-file parser)
                                                  sum (min num max-sequences-per-file)))
                                  (loop for parser in parsers
                                        sum (loop for num in (num-sequences-per-file parser)
                                                  sum num))))
         (_ (format t "estimated sequences: ~a~%" estimated-sequences))
         (number-of-updates 200)
         (update-increment (floor estimated-sequences number-of-updates))
         (next-update update-increment)
         (tracker (make-instance 'job-tracker))
         (number-of-resources (+ 5 (* number-of-workers 10)))
         (resources (cons (loop for index to (* 2 number-of-resources)
                                collect (make-sequence-align))
                          nil))
         (queue (lparallel.queue:make-queue :fixed-capacity number-of-resources))
         (workers (loop for index below number-of-workers
                        for worker = (make-instance 'worker
                                                    :worker-index index
                                                    :tracker tracker
                                                    :resource-stack resources)
                        collect worker)))
    ;; Reader process
    (vformat verbose t "About to start reader~%")
    (loop for parser in parsers
          for threads = (loop for worker in workers
                              for thread = (bt:make-thread
                                            (lambda ()
                                              (vformat *verbose* t "reader About to start thread~%")
                                              (wait-for-sequence *worker*
                                                                 *queue*
                                                                 :verbose verbose))
                                            :initial-bindings (list (cons '*worker* worker)
                                                                    (cons '*queue* queue)
                                                                    (cons '*verbose* verbose)))
                              collect thread)
          for results = (make-hash-table :test 'equal :thread-safe t)
          for results-lock = (bt:make-lock "results-lock")
          for parser-sequence-count = 0
          for parser-name = (name parser)
          for files = (files parser)
          do (unwind-protect
                  (let* ((pass-fail  (when pass-fail-files
                                       (let* ((pass-fail-dir (make-pathname :name nil :type nil :defaults (first files)))
                                              (pass-file-name (make-pathname
                                                               :name (format nil "~a-pass" parser-name)
                                                               :type "fastq" :defaults pass-fail-dir))
                                              (fail-file-name (make-pathname
                                                               :name (format nil "~a-fail" parser-name)
                                                               :type "fastq" :defaults pass-fail-dir)))
                                         (make-instance 'pass-fail
                                                        :pass-file (sa:make-seq-file-out (namestring pass-file-name))
                                                        :pass-file-lock (bt:make-lock "pass-file-lock")
                                                        :fail-file (sa:make-seq-file-out (namestring fail-file-name))
                                                        :fail-file-lock (bt:make-lock "fail-file-lock"))))))
                    (unwind-protect
                         (loop named files-loop
                               for file in files
                               for seq-file = (seqan:make-seq-file-in (namestring file))
                               for file-sequence-count = 0
                               do (unwind-protect
                                       (loop named read-loop
                                             for data = (progn
                                                          (vformat verbose t "reader: About to pop resource~%")
                                                          (mp:atomic-pop (car resources)))
                                             do (progn
                                                  (vformat verbose t "reader: About to read a record into data: ~a~%" data)
                                                  (when (sa:at-end seq-file)
                                                    (return-from read-loop nil))
                                                  (sa:read-record (metadata data)
                                                                  (seq data)
                                                                  seq-file)
                                                  (setf (pass-fail data) pass-fail
                                                        (results-lock data) results-lock
                                                        (results data) results)
                                                  (incf sequence-count)
                                                  (incf file-sequence-count)
                                                  (when (and progress-callback (>= sequence-count next-update))
                                                    (incf next-update update-increment)
                                                    (multiple-value-bind (progress msg)
                                                        (evaluate-progress start-time
                                                                           (mp:atomic (slot-value tracker 'sequence-count))
                                                                           (mp:atomic (slot-value tracker 'good-sequence-count))
                                                                           estimated-sequences
                                                                           file-sequence-count
                                                                           file)
                                                      (funcall progress-callback progress msg)))
                                                  (when (and max-sequences-per-file (>= file-sequence-count max-sequences-per-file))
                                                    (return-from read-loop nil))
                                                  (vformat verbose t "reader: About to push-queue ~a~%" (sa:to-string (seq data)))
                                                  (if (> number-of-workers 0)
                                                      (lparallel.queue:push-queue data queue)
                                                      (mp:atomic-push data (car resources))))) ; with no workers just push resource back
                                    (sa:close seq-file)))
                      (when pass-fail
                        (sa:close (pass-file pass-fail))
                        (sa:close (fail-file pass-fail)))))
               ;; Shutdown the workers
               (progn
                 (vformat verbose t "reader: sending quit message to ~a threads~%" (length threads))
                 (loop for thread in threads
                       do (unless testing
                            (lparallel.queue:push-queue :quit queue)))
                 ;; Wait for them to shutdown
                 (vformat verbose t "reader: waiting for ~a threads to quit~%" (length threads))
                 (loop for thread in threads
                       do (bt:join-thread thread)
                       do (vformat verbose t "Thread ~a shutdown~%" thread))))
          do ;; Accumulate the results
             (progn
               (vformat verbose t "reader: Accumulating worker results~%")
               (format t "Read ~a sequences~%" sequence-count)
               (format t "Results: ~s~%" results)
               ;; in result
               (let* ((filtering (make-instance 'filtering
                                                :name (name parser)
                                                :file-names (files parser)
                                                :sequences results
                                                :start-time start-time
                                                :current-time (get-internal-real-time)
                                                :total-sequence-count parser-sequence-count)))
                 (format t "Writing results to ~a~%" (output-file-name parser))
                 (save-filtering filtering (output-file-name parser) :overwrite (overwrite parser)))))
    (format t "Wrote results to: ~a~%" (mapcar 'output-file-name parsers))
    ))

(defun make-job-tracker-vector (number-of-jobs)
  (let ((trackers (make-array number-of-jobs)))
    (loop for index below number-of-jobs
          do (setf (elt trackers index) (make-instance 'job-tracker)))
    trackers))

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
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

|#
