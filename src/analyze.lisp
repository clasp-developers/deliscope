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
  ((id :initarg :id :accessor id)
   (metadata :initarg :metadata :accessor metadata)
   (cl-string :initform (make-array 170 :element-type 'base-char
                                        :adjustable t
                                        :fill-pointer 0
                                        :initial-element #\*)
              :accessor cl-string
              :documentation "Preallocated CL string to spill sequence into for string ops")
   (seq :initarg :seq :accessor seq)
   (align :initarg :align :accessor align)
   (row0& :initarg :row0& :reader row0&)
   (row1& :initarg :row1& :reader row1&)
   (good-result-key :initarg :good-result-key :accessor good-result-key)))

(defparameter *seqan-make-lock* (bt:make-lock 'seqan-make-lock))
(defun make-sequence-align (&optional (id 0))
  "seqan doesn't appear to have thread safe allocators - so we use a lock
when allocating align and string objects"
  (unwind-protect
       (progn
         (bt:acquire-lock *seqan-make-lock*)
         (let ((align (seqan:make-align :dna5q-string)))
           (seqan:resize (sa:rows& align) 2)
           (make-instance 'sequence-align
                          :id id
                          :metadata (seqan:make-string :char-string)
                          :seq (seqan:make-string :dna5q-string)
                          :align align
                          :row0& (seqan:row& align 0)
                          :row1& (seqan:row& align 1))))
    (bt:release-lock *seqan-make-lock*)))

(defmethod print-object ((object sequence-align) stream)
  (print-unreadable-object (object stream :type t)
    (format stream "~a" (id object))))

(defclass parser ()
  ((name :initarg :name :accessor name)
   (files :initarg :files :accessor files)
   (num-sequences-per-file :initarg :num-sequences-per-file :accessor num-sequences-per-file)
   (output-file-name :initarg :output-file-name :accessor output-file-name)
   (overwrite :initarg :overwrite :accessor overwrite)
   ))

(defclass review ()
  ((name :initarg :name :accessor name)
   (file-names :initarg :file-names :accessor file-names)
   (sequences :initform (make-hash-table :test 'equal) :initarg :sequences :accessor sequences)
   (start-time :initarg :start-time :accessor start-time)
   (current-time :initarg :current-time :accessor current-time)
   (total-read-count :initform 0 :initarg :total-read-count :accessor total-read-count)
   (good-read-count :initform 0 :initarg :good-read-count :accessor good-read-count)))

(defmethod summarize (review)
  (format t "Files:~%~s~%" (file-names review))
  (format t "Total read: ~a~%" (total-read-count review))
  (format t "Good read:  ~a~%" (good-read-count review))
  (format t "Sample sequences:~%")
  (show-hash-table (sequences review))
  (values))

(defclass worker ()
  ((resource-stack :initarg :resource-stack :reader resource-stack)
   (tracker :initarg :tracker :reader tracker)
   (worker-index :initarg :worker-index :reader worker-index)
   (worker-log :initarg :worker-log :reader worker-log)
   (thread :initform nil :accessor thread)))


(defclass batch ()
  ((pass-fail-data :initarg :pass-fail-data :accessor pass-fail-data)
   (results-lock :initarg :results-lock :accessor results-lock)
   (results-data :initarg :results-data :accessor results-data)
   (batch-list :initarg :batch-list :accessor batch-list)))

(defclass job-tracker ()
  ((total-read-count :initform 0 :accessor total-read-count)
   (good-read-count :initform 0 :accessor good-read-count)
   (parser-total-read-count :initform 0 :accessor parser-total-read-count)
   (parser-good-read-count :initform 0 :accessor parser-good-read-count)))

(defun zero-parser-counts (job-tracker)
  (setf (parser-total-read-count job-tracker) 0
        (parser-good-read-count job-tracker) 0))

(defmethod print-object ((obj job-tracker) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream ":total-read-count ~a :good-read-count: ~a :parser-total-read-count ~a :parser-good-read-count: ~a"
            (total-read-count obj)
            (good-read-count obj)
            (parser-total-read-count obj)
            (parser-good-read-count obj)
            )))

(defun sum-trackers (sum-tracker trackers)
  (loop for index below (length trackers)
        for tracker = (elt trackers index)
        sum (total-read-count tracker) into total-read-count-sum
        sum (good-read-count tracker) into good-read-count-sum
        finally (setf (total-read-count sum-tracker) total-read-count-sum
                      (good-read-count sum-tracker) good-read-count-sum)))

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
    (format stream "~a ~a :low-quality-count ~a ~s"
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

(defmacro wformat (worker fmt &rest args)
  `(when (worker-log ,worker)
     (format (worker-log ,worker) "w~a: " (worker-index worker))
     (format (worker-log ,worker) ,fmt ,@args)))

(defmacro wformat* (worker fmt &rest args)
  `(when (worker-log ,worker)
     (format (worker-log ,worker) ,fmt ,@args)))

(defun analyze-column-using-align (worker sequence-align column-codes start adjust low-quality-count quality-string)
  (let* ((abs-start (+ start (start column-codes) adjust))
         (best-score -999999)
         (best-row-code nil))
    (loop for row-code in (codes column-codes)
          for name = (del-name row-code)
          for code = (code row-code)
          for lower-diag = (- abs-start 1)
          for upper-diag = (+ lower-diag 1) ;;+ start end-offset )
          for score = (progn
                        (sa:assign-source (row0& sequence-align) (seq sequence-align))
                        (sa:assign-source (row1& sequence-align) code)
                        (sa:global-alignment (align sequence-align) *simple-score* *align-config*
                                                          lower-diag
                                                          upper-diag ))
          do (wformat worker "Score: ~a~%" score)
             (wformat* worker "~a~%" (sa:to-string (align sequence-align)))
          do (when (> score best-score)
               (setf best-score score
                     best-row-code row-code)))
    (make-instance 'match
                   :score best-score
                   :low-quality-count low-quality-count
                   :quality-string quality-string
                   :row-code best-row-code)))

(defun analyze-column-fast (worker sequence-align column-codes start adjust)
  (seqan:assign (cl-string sequence-align) (seq sequence-align))
  (let* ((seq-string (cl-string sequence-align))
         (abs-start (+ start (start column-codes) adjust))
         (abs-end (+ start (end column-codes) adjust)))
    (wformat worker "abs-start: ~d  abs-end ~d Subsequence: ~a~%" abs-start abs-end (subseq (cl-string sequence-align) (1- abs-start) (1- abs-end)))
    (loop for row-code in (codes column-codes)
          for code-seq = (sequence row-code)
          when (string= seq-string code-seq :start1 (1- abs-start) :end1 (1- abs-end))
            do (progn
                 (wformat worker "Fast matched with row-code ~a~%" row-code)
                 (return-from analyze-column-fast (values t row-code))))
    (values nil)))

(defun analyze-column (worker sequence-align sequence-length column-codes start adjust)
  (let* ((abs-start (+ start (start column-codes) adjust))
         (abs-end (+ start (end column-codes) adjust)))
    (wformat worker "abs-end = ~a   sequence-length ~a~%" abs-end sequence-length)
    (cond
      ((<= abs-end sequence-length)
        (multiple-value-bind (low-quality-count qual-string)
            (sa:count-quality-value-less-than (seq sequence-align) *minimum-phred-quality* abs-start abs-end)
          (multiple-value-bind (fast-worked row-code)
              (analyze-column-fast worker sequence-align column-codes start adjust)
            (when fast-worked
              (return-from analyze-column (make-instance 'match :score 0
                                                                :low-quality-count low-quality-count
                                                                :row-code row-code
                                                                :quality-string qual-string))))
          (analyze-column-using-align worker sequence-align column-codes start adjust low-quality-count qual-string)))
      (t 
       (wformat worker "Missing code for ~a - ~a sequence stops at ~a~%" abs-start abs-end sequence-length)
       :incomplete-code))))

(defun analyze-forward-primer (worker sequence-align)
  (sa:assign-source (row0& sequence-align) (seq sequence-align))
  (sa:assign-source (row1& sequence-align) (forward-primer *all-codes*))
  (let ((fwd-score (sa:global-alignment (align sequence-align) *simple-score* *align-config*))
        (start (1+ (sa:to-view-position (row1& sequence-align) (sa:length (forward-primer *all-codes*))))))
    (wformat worker "Alignment: ~%~a~%" (sa:to-string (align sequence-align)))
    (values start fwd-score)))

(defun analyze-sequence (worker sequence-align)
  (let ((sequence (seq sequence-align)))
    (wformat worker "analyze-sequence for ~a start with forward-primer~%" sequence-align)
    (multiple-value-bind (start fwd-score)
        (analyze-forward-primer worker sequence-align)
      (wformat worker "(length fwd-primer): ~a~%" (sa:length (forward-primer *all-codes*)))
      (wformat worker "fwd-score: ~a (0 best) 3' end: ~a~%" fwd-score start)
      (let ((length-sequence (sa:length sequence))
            (length-codes (+ start (end *all-codes*))))
        (wformat worker "length-sequence -> ~a   length-codes -> ~a~%" length-sequence length-codes)
        (let ((result (loop for column-codes in (column-codes *all-codes*)
                            for start-offset = (start column-codes)
                            for end-offset = (end column-codes)
                            for best-digit = (analyze-column worker sequence-align length-sequence column-codes start 0)
                            do (wformat worker "offset: ~a best-digits: ~a~%" start-offset best-digit)
                            when (not (eq best-digit :incomplete-code))
                              collect best-digit)))
          result)))))

(defvar *review* (make-hash-table :test 'equal))
(defvar *bead-specific (make-hash-table :test 'equal))


(defun accept-result (result)
  (and result
       (every (lambda (x) (>= (score x) -2)) result)
       (every (lambda (x) (< (low-quality-count x) *low-quality-max*)) result)))

(defun evaluate-progress (start-time count good-read-count total-count file-read-count file)
  (if (= count 0)
      (let ((msg (format nil "Read: ~a from ~a~%Parsed ~a of ~a~%Estimating remaining time...~%" file-read-count file count total-count)))
        (values 0 msg))
      (let* ((current-time (get-internal-real-time))
             (elapsed-seconds (float (/ (- current-time start-time) internal-time-units-per-second)))
             (elapsed-minutes (/ elapsed-seconds 60.0))
             (elapsed-hours (/ elapsed-minutes 60.0))
             (total-seconds (if (= count 0) 99999999.0 (* (/ total-count count) elapsed-seconds)))
             (remaining-seconds (- total-seconds elapsed-seconds)) ; overestimate
             (remaining-minutes (/ remaining-seconds 60.0))
             (remaining-hours (/ remaining-minutes 60.0)))
        (multiple-value-bind (safe-elapsed-hours safe-elapsed-minutes)
            (if (< elapsed-hours 0.0)
                (values 0 0)
                (floor elapsed-hours))
          (multiple-value-bind (hours minutes)
              (if (< remaining-hours 0.0)
                  (values 0 0)
                  (floor remaining-hours))
            (labels ((time-description (hours minutes)
                       (if (zerop hours)
                           (format nil "~4,1f minute~:p" minutes)
                           (format nil "~a hour~:p ~4,1f minute~:p" hours minutes))))
              (let ((msg (format nil "Read ~a sequences from ~a~%Parsed ~a of ~a total reads / ~a good reads~%Elapsed time:~a   Estimated remaining time:~a~%"
                                 file-read-count file
                                 count total-count
                                 good-read-count
                                 (time-description safe-elapsed-hours (* 60.0 safe-elapsed-minutes))
                                 (time-description hours (* 60.0 minutes))
                                 ))
                    (progress (floor (float (* 100.0 (/ count total-count))))))
                (values progress msg))))))))


(defun wait-for-sequence (worker queue &key verbose)
  (handler-case
      (loop named work-loop
            for batch = (let ((bbb (progn
                                     (wformat worker "About to pop queue queue-count: ~a~%" (lparallel.queue:queue-count queue))
                                     (lparallel.queue:pop-queue queue))))
                          (when (eq bbb :quit)
                            (wformat worker "exiting~%")
                            (return-from work-loop nil))
                          (unless bbb
                            (wformat worker "got nil as bbb~%"))
                          bbb)
            for total-read-count = 0
            for good-read-count = 0
            do (wformat worker "got batch: ~a batch-list-len ~a~%" batch (length (batch-list batch)))
            do (loop named batch-loop
                     for data in (batch-list batch)
                     for __1 = (wformat worker "data: ~a~%" data)
                     for metadata = (metadata data)
                     for seq = (seq data)
                     for result = (analyze-sequence worker data)
                     for _1 = (wformat worker "result = ~a~%" result)
                     for good-result-or-nil = (if (accept-result result)
                                                  (progn
                                                    (wformat worker "Accepted result~%")
                                                    (incf good-read-count)
                                                    (mapcar (lambda (digit) (when digit (del-name (row-code digit)))) result))
                                                  (progn
                                                    (wformat worker "Discarded result~%")

                                                    nil))
                     do (incf total-read-count)
                     do (setf (good-result-key data) good-result-or-nil))
            do (wformat worker "Inc sequence count~%")
            do (mp:atomic-incf (slot-value (tracker worker) 'total-read-count) total-read-count)
            do (mp:atomic-incf (slot-value (tracker worker) 'parser-total-read-count) total-read-count)
            do (when (> good-read-count 0)
                 (mp:atomic-incf (slot-value (tracker worker) 'good-read-count) good-read-count)
                 (mp:atomic-incf (slot-value (tracker worker) 'parser-good-read-count) good-read-count)
                 (bt:with-lock-held ((results-lock batch))
                   (loop named record-results
                         for data in (batch-list batch)
                         for good-result-key = (good-result-key data)
                         when good-result-key
                           do (incf (gethash good-result-key (results-data batch) 0))))
                 (when (pass-fail-data batch)
                   (bt:with-lock-held ((pass-file-lock (pass-fail-data batch)))
                     (sa:write-record (pass-file (pass-fail-data batch))))))
            do (wformat worker "maybe write fail~%")
            do (when (/= good-read-count total-read-count)
                 (when (pass-fail-data batch)
                   (bt:with-lock-held ((fail-file-lock (pass-fail-data batch)))
                     (sa:write-record (fail-file (pass-fail-data batch))))))
            do (wformat worker "recycle~%")
            do (loop named recycle-loop
                     for data in (batch-list batch)
                     do (wformat worker "Recycling ~a~%" data)
                     do (mp:atomic-push data (car (resource-stack worker))) ; recycle the data 
                     ))
    (error (c)
      (if verbose
          (wformat worker "Error wait-for-sequence:  ~a~%" c)
          (format t "Error wait-for-sequence: ~a~%" c)))))


(defun maybe-make-pass-fail (pass-fail-files files parser-name)
  (when pass-fail-files
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
                     :fail-file-lock (bt:make-lock "fail-file-lock")))))

(defun read-batch (batch-size seq-file resources verbose)
  (declare (ignorable verbose))
  (loop named read-batch
        with batch = nil
        for data-index below batch-size
        collect (if (sa:at-end seq-file)
                    (return-from read-batch batch)
                    (let* ((data (mp:atomic-pop (car resources))))
                      #+(or)(when verbose (format t "Reader about to read into data: ~a~%" data))
                      (sa:read-record (metadata data)
                                      (seq data)
                                      seq-file)
                      data))))

(defun update-progress (progress-callback start-time tracker total-read-count estimated-sequences file-read-count file)
  (multiple-value-bind (progress msg)
      (evaluate-progress start-time
                         total-read-count #+(or)(mp:atomic (slot-value tracker 'total-read-count))
                         (mp:atomic (slot-value tracker 'good-read-count))
                         estimated-sequences
                         file-read-count
                         file)
    (funcall progress-callback progress msg)))

(defun shutdown-and-wait-for-workers (threads testing queue verbose)
  (when verbose (format t "reader: sending quit message to ~a threads~%" (length threads)))
  (loop for thread in threads
        do (unless testing
             (lparallel.queue:push-queue :quit queue)))
  ;; Wait for them to shutdown
  (when verbose (format t "reader: waiting for ~a threads to quit~%" (length threads)))
  (loop for thread in threads
        do (bt:join-thread thread)
        do (format t "Thread ~a shutdown~%" thread)))

(defparameter *worker* nil)
(defparameter *queue* nil)
(defparameter *verbose* nil)
(defun analyze-parsers-using-workers (number-of-workers parsers
                                      &key max-sequences-per-file progress-callback pass-fail-files verbose testing)
  (format t "Entered analyze-files-using-workers ~a~%" parsers)
  ;; Create resources and workers
  (format t "max-sequences-per-file: ~a~%" max-sequences-per-file)
  (let* ((total-read-count 0)
         (start-time (get-internal-real-time))
         (estimated-sequences (if max-sequences-per-file
                                  (loop for parser in parsers
                                        sum (loop for num in (num-sequences-per-file parser)
                                                  sum (min num max-sequences-per-file)))
                                  (loop for parser in parsers
                                        sum (loop for num in (num-sequences-per-file parser)
                                                  sum num))))
         (_1 (format t "estimated sequences: ~a~%" estimated-sequences))
         (number-of-updates 200)
         (batch-size 100)
         (queue-size (+ number-of-workers 100))
         (update-increment (floor estimated-sequences number-of-updates))
         (next-update update-increment)
         (tracker (make-instance 'job-tracker))
         (number-of-resources (+ 10 (* queue-size batch-size)))
         (resources (cons (loop for index to (* 2 number-of-resources)
                                collect (make-sequence-align index))
                          nil))
         (queue (lparallel.queue:make-queue :fixed-capacity queue-size))
         (workers (loop for index below number-of-workers
                        for worker = (make-instance 'worker
                                                    :worker-index index
                                                    :worker-log (if verbose
                                                                    (open (format nil "worker~a.log" index) :direction :output :if-exists :supersede)
                                                                    nil)
                                                    :tracker tracker
                                                    :resource-stack resources)
                        collect worker)))
    (declare (ignore _1))
    ;; Reader process
    (format t "Batch size: ~a~%" batch-size)
    (loop for parser in parsers
          for threads = (loop for worker in workers
                              for thread = (progn
                                             (zero-parser-counts (tracker worker))
                                             (bt:make-thread
                                              (lambda ()
                                                (wformat worker "reader About to start thread~%")
                                                (wait-for-sequence *worker*
                                                                   *queue*
                                                                   :verbose verbose))
                                              :initial-bindings (list (cons '*worker* worker)
                                                                      (cons '*queue* queue)
                                                                      (cons '*verbose* verbose))))
                              collect thread)
          for results = (make-hash-table :test 'equal :thread-safe t)
          for results-lock = (bt:make-lock "results-lock")
          for parser-name = (name parser)
          for files = (files parser)
          do (unwind-protect
                  (let ((pass-fail (maybe-make-pass-fail pass-fail-files files parser-name)))
                    (unwind-protect
                         (loop named files-loop
                               for file in files
                               for seq-file = (seqan:make-seq-file-in (namestring file))
                               for file-read-count = 0
                               do (unwind-protect
                                       (loop named read-loop
                                             do (when (sa:at-end seq-file)
                                                  (return-from read-loop nil))
                                             do (let ((data (read-batch batch-size seq-file resources verbose)))
                                                  (incf total-read-count (length data))
                                                  (incf file-read-count (length data))
                                                  #+(or)(format t "Reader: About to push-queue queue count: ~a~%" (lparallel.queue:queue-count queue))
                                                  (lparallel.queue:push-queue (make-instance 'batch
                                                                                             :pass-fail-data pass-fail
                                                                                             :results-lock (bt:make-lock "results-lock")
                                                                                             :results-data results
                                                                                             :batch-list data)
                                                                              queue)
                                                  (when (and progress-callback (>= total-read-count next-update))
                                                    (incf next-update update-increment)
                                                    (update-progress progress-callback start-time tracker total-read-count estimated-sequences file-read-count file))
                                                  (when (and max-sequences-per-file (>= file-read-count max-sequences-per-file))
                                                    (return-from read-loop nil))))
                                    (sa:close seq-file)))
                      (when pass-fail
                        (sa:close (pass-file pass-fail))
                        (sa:close (fail-file pass-fail)))))
               ;; Shutdown the workers
               (shutdown-and-wait-for-workers threads testing queue verbose))
          do ;; Accumulate the results
             (progn
               (format t "reader: Accumulating worker results~%")
               (format t "Read ~a sequences~%" total-read-count)
               (format t "Results: ~s~%" results)
               (format t "Tracker: ~s~%" tracker)
               ;; in result
               (let* ((review (make-instance 'review
                                                :name (name parser)
                                                :file-names (files parser)
                                                :sequences results
                                                :start-time start-time
                                                :current-time (get-internal-real-time)
                                                :total-read-count (parser-total-read-count tracker)
                                                :good-read-count (parser-good-read-count tracker))))
                 (format t "Writing results to ~a~%" (output-file-name parser))
                 (save-review review (output-file-name parser) :overwrite (overwrite parser)))))
    (format t "Wrote results to: ~a~%" (mapcar 'output-file-name parsers))
    ))

(defun make-job-tracker-vector (number-of-jobs)
  (let ((trackers (make-array number-of-jobs)))
    (loop for index below number-of-jobs
          do (setf (elt trackers index) (make-instance 'job-tracker)))
    trackers))

(defun sequence-survey (seq review)
  (let ((ht (make-hash-table :test 'equal)))
    (maphash (lambda (k v)
               (let ((seq-key (subseq k 1 4)))
                 (push (cons k v) (gethash seq-key ht))))
             review)
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

(defun unique-library-codes (sequences)
  "Gather up all of the unique library codes and count how many times they appear"
  (let ((library-codes (make-hash-table :test 'equal)))
    (maphash (lambda (key value)
               (let ((library-code (subseq key 4 5)))
                 (incf (gethash (first library-code) library-codes 0) value)))
             sequences)
    library-codes))

(defun bead-codes (review &optional (value-cutoff 150))
  "Gather the sequences and count how many times they appear with different bead codes"
  (let ((bead-codes (make-hash-table :test 'equal)))
    (maphash (lambda (key value)
               (let ((bead-code (subseq key 0 2))
                     (seq-code (subseq key 2 8)))
                 (when (>= value value-cutoff)
                   (push (list bead-code value) (gethash seq-code bead-codes)))))
             (sequences review))
    bead-codes))

(defun bin-codes (review)
  (let ((digit-counts (make-hash-table :test 'equal)))
    (maphash (lambda (key value)
               (declare (ignore value))
               (loop for digit-index below 8
                     for code = (copy-seq (elt key digit-index))
                     do (setf (elt code 1) #\_)
                     do (incf (gethash code digit-counts 0))))
             (sequences review))
    digit-counts))

(defun safe-ratio (top bot)
  (if (= bot 0)
      -9999.9
      (/ (float top) (float bot))))

(defun describe-review (review)
  (format t "Loaded ~a sequences for ~a~%" (hash-table-count (sequences review)) (name review))
  (format t "   There were ~a total reads and ~a were accepted = ~5,2f%~%"
          (total-read-count review)
          (good-read-count review)
          (* 100.0 (safe-ratio (good-read-count review) (total-read-count review))))
  (let ((raw-library-codes (make-hash-table :test 'equal))
        (short-keys (make-hash-table :test 'eql)))
    (maphash (lambda (key value)
               (declare (ignore value))
               (when (< (length key) 10)
                 (incf (gethash (length key) short-keys 0)))
               (when (= (length key) 10)
                 (let ((library-key (subseq key 8 10)))
                   (incf (gethash library-key raw-library-codes 0)))))
             (sequences review))
    (when (> (hash-table-count short-keys) 0)
      (format t "   Short keys:~%")
      (maphash (lambda (k v)
                 (format t "     length ~a num ~a~%" k v))
               short-keys))
    (format t "   Raw library-codes:~%")
    (maphash (lambda (k v)
               (format t "     ~s ~a~%" k v))
             raw-library-codes))
  )

(defun save-review (review file &key overwrite)
  "Save the results to a file. Pass :force t if you want to force overwriting an existing file"
  (when (and (not overwrite) (probe-file file))
    (error "The file ~a already exists~%" file))
  (with-open-file (fout file :direction :output :if-exists (if overwrite :supersede))
    (let ((*print-readably* t)
          (*print-pretty* nil))
      (prin1 (name review) fout)
      (terpri fout)
      (prin1 (file-names review) fout)
      (terpri fout)
      (prin1 (total-read-count review) fout)
      (terpri fout)
      (prin1 (good-read-count review) fout)
      (terpri fout)
      (prin1 (hash-table-count (sequences review)) fout)
      (terpri fout)
      (maphash (lambda (k v)
                 (format fout "~s ~s~%" k v))
               (sequences review))))
  (format t "Saved to ~a~%" file))

(defun load-review (file &key (verbose t))
  (with-open-file (fin (merge-pathnames file) :direction :input)
    (let ((review (make-instance 'review)))
      (setf (name review) (read fin))
      (setf (file-names review) (read fin))
      (setf (total-read-count review) (read fin))
      (setf (good-read-count review) (read fin))
      (let* ((num-pairs (read fin)))
        (loop for index below num-pairs
              for key = (read fin)
              for val = (read fin)
              do (setf (gethash key (sequences review)) val))
        (when verbose (describe-review review))
        review))))

(defun fix-library-code (review library-x library-y)
  "For each key fix up the library code.
If library-x is missing then append library-x and library-y.
if library-x matches the key and library-y is missing then append library-y.
If both library codes are present kep the key"
  (let ((new-seqs (make-hash-table :test 'equal)))
    (maphash (lambda (key value)
               (cond
               ((= (length key) 10)
                (setf (gethash key new-seqs) value))
               ((and (= (length key) 9)
                     (string= (elt key 8) library-x))
                (setf (gethash (append key (list library-y)) new-seqs) value))
               ((= (length key) 8)
                (setf (gethash (append key (list library-x library-y)) new-seqs) value))
               (t );discard
               ))
             (sequences review))
    (setf (sequences review) new-seqs))
  review)

(defun join-codons (input)
  (let ((review (make-hash-table :test 'equal)))
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
                 (setf (gethash compact-sequence-key review) val)))
             input)
    review))

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
               (when (and (string= (elt key 4) *library-code*)
                          (> value threshold))
                 (setf (gethash key result-sequences) value)))
             input)
    result-sequences))

(defun bead-specific-counts (input-review &key (key-positions '(1 2 3)) library-code)
  (let ((review (make-hash-table :test 'equal)))
    (maphash (lambda (key val)
               (declare (ignore val))
               (let ((library-key (elt key 4))
                     (sequence-key (loop for pos in key-positions
                                         collect (elt key pos)))) 
                 (if library-code
                     (when (string= library-code library-key)
                       (incf (gethash sequence-key review 0)))
                     (incf (gethash sequence-key review 0)))))
             input-review)
    review))

(defclass analysis (review)
  ((joined-sequences :initarg :jointed-sequences :accessor joined-sequences)
   (minimum-counts :initform *minimum-counts* :accessor minimum-counts)
   (filtered-sequences :accessor filtered-sequences)
   (library-code :initform *library-code* :accessor library-code)
   (keys :initarg :keys :accessor keys)
   (redundancy :accessor redundancy)
   ))

(defclass simple-hit ()
  ((sort-day :initarg :sort-day :accessor sort-day)
   (num-reads :initarg :num-reads :accessor num-reads)))

(defmethod print-object ((obj simple-hit) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream ":sort-day ~s :num-reads ~a"
            (sort-day obj)
            (num-reads obj))))

(defclass hit (simple-hit)
  ((bead-code :initform nil :initarg :bead-code :accessor bead-code)))

(defmethod print-object ((obj hit) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream ":sort-day ~s :num-reads ~a :bead-code ~s"
            (sort-day obj)
            (num-reads obj)
            (bead-code obj))))

(defclass conclusion ()
  ((results :initarg :results :accessor results)))

(defun process-results (review &key (keys '(1 2 3)))
  (change-class review 'analysis :keys keys)
  (setf (joined-sequences review) (join-codons (sequences review)))
  (setf (filtered-sequences review) (filter-sequences (joined-sequences review) (minimum-counts review)))
  (setf (redundancy review) (bead-specific-counts (filtered-sequences review)
                                                     :key-positions keys
                                                     :library-code (library-code review)))
  review)


(defun filtered-library-codes (analysis)
  (let ((lib-codes (unique-library-codes (filtered-sequences analysis))))
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
(defun sort-hits (&optional review)
  (let ((unsorted))
    (maphash (lambda (key value)
               (push (cons key value) unsorted))
             review)
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



(defclass parsers ()
  ((parsers :initform (make-hash-table :test 'eq :thread-safe t)
            :initarg :parsers
            :accessor parsers))
  )

(defun estimate-sequences-in-file (file-name &key num verbose)
  (declare (ignore verbose))
  (unless (probe-file file-name)
    (error "Could not find file: ~a" file-name))
  (if num
      num
      (multiple-value-bind (num-lines file-pos file-size)
          (core:count-lines-in-file (namestring file-name) 1000000)
          (let ((estimate-lines (floor (* (/ file-size file-pos) num-lines))))
            (floor (/ estimate-lines 4))))))


(defvar *data-directory* nil)
(defun set-data-directory (directory)
  (setf directory (uiop:ensure-directory-pathname (pathname directory)))
  (unless (probe-file directory)
    (error "The directory ~a does not exist" directory))
  (format t "data directory set to ~a~%" directory)
  (setf *data-directory* directory))

(defun create-sequence-parser (name &rest args)
  (let ((output-file-name (make-pathname :name (format nil "~a-results" (string-downcase (string name)))
                                         :type "dat"))
        (overwrite t))
    (when (eq (car args) :output)
      (pop args)
      (setf output-file-name (pop args)))
    (let* ((files args)
           (absolute-files (mapcar (lambda (file)
                                     (if *data-directory*
                                         (merge-pathnames file *data-directory*)
                                         file))
                                   files)))
      (format t "Scanning the files to count the number of sequences.~%")
      (format t "This may take a few minutes - hit ii to interrupt.~%")
      (finish-output)
      (let* ((num-seq-per-file (mapcar (lambda (file)
                                         (estimate-sequences-in-file file :verbose nil))
                                       absolute-files))
             (parser (make-instance 'parser :name name
                                            :files absolute-files
                                            :num-sequences-per-file num-seq-per-file
                                            :output-file-name output-file-name
                                            :overwrite overwrite)))
        (format t "Number of sequences: ~a~%" (apply '+ num-seq-per-file))
        (finish-output)
        parser))))





(defun parse-tty (parsers &key (limit nil))
  "Usage: (parse-tty (list <parsers>...) [:limit <integer>] )

Parse the sequences specified by each parser.
You can limit the number sequences loaded from each file within each parser by adding :limit <num>.
"
  (declare (ignore action instance))
  (format t "Number of threads: ~a~%" (core:num-logical-processors))
  (analyze-parsers-using-workers
   (core:num-logical-processors)
   parsers
   :max-sequences-per-file limit
   :progress-callback (let ((last-val -1))
                        (lambda (val msg &key done)
                          (if done
                              (format t "Done.~%")
                              (progn
                                (when (> val last-val)
                                  (setf last-val val)
                                  (progn
                                    (format t "Progress: ~a~%" val)
                                    (format t "~a~%" msg)
                                    (finish-output)
                                    )
                                  )
                                )))))
  (values))

(defun draw-conclusion-with-bead-codes (analyses)
  (let ((result-ht (make-hash-table :test 'equalp)))
    (loop for analysis in analyses
          do (maphash (lambda (key value)
                        (let* ((bead-code (elt key 0))
                               (agg-key (subseq key 1 4))
                               (agg-value (make-instance 'hit
                                                         :sort-day (name analysis)
                                                         :num-reads value
                                                         :bead-code bead-code)))
                          (let ((rht (alexandria:ensure-gethash agg-key result-ht (make-hash-table :test 'equalp))))
                            (setf (gethash (list (name analysis) bead-code) rht) agg-value))))
                      (filtered-sequences analysis)))
    (make-instance 'conclusion
                   :results result-ht)))


(defun redraw-conclusion-ignoring-bead-codes (conclusion)
  (let ((results-ht (make-hash-table :test 'equalp)))
    (maphash (lambda (key hits)
               (let ((sort-ht (make-hash-table)))
                 (maphash (lambda (key hit)
                            (declare (ignore key))
                            (let ((simple-hit (gethash (sort-day hit) sort-ht)))
                              (if simple-hit
                                  (progn
                                    (incf (num-reads hit) (num-reads simple-hit))
                                    (push (bead-code hit) (bead-code simple-hit)))
                                  (let ((new-simple-hit (make-instance 'hit
                                                                       :sort-day (sort-day hit)
                                                                       :num-reads (num-reads hit)
                                                                       :bead-code (list (bead-code hit)))))
                                    (setf (gethash (sort-day hit) sort-ht) new-simple-hit)))))
                          hits)
                 (setf (gethash key results-ht) sort-ht)))
             (results conclusion))
    (make-instance 'conclusion
                   :results results-ht)))

(defun write-conclusion-to-stream (conclusion stream)
  (let ((unsorted-results nil))
    (maphash (lambda (key value)
               (push (cons key value) unsorted-results))
             (results conclusion))
    (let ((results (sort unsorted-results #'> :key (lambda (v) (hash-table-count (cdr v))))))
      (loop for result in results
            for key = (car result)
            for values = (cdr result)
            for index from 1
            when (> (hash-table-count values) 1)
              do (progn
                   (format stream "#~a ~a seq: ~s~%" index (hash-table-count values) key)
                   (maphash (lambda (key hit)
                              (declare (ignore key))
                              (if (bead-code hit)
                                  (format stream "    sort-day(~a) bc(~s) reads(~8d)~%"
                                          (sort-day hit) (bead-code hit) (num-reads hit))
                                  (format stream "    sort-day(~a) reads(~8d)~%"
                                          (sort-day hit) (num-reads hit))
                                  ))
                            values))))))

(defun describe-conclusion (conclusion)
  (write-conclusion-to-stream conclusion t))

(defun write-conclusion-to-file (conclusion filename)
  (with-open-file (fout filename :direction :output)
    (write-conclusion-to-stream conclusion fout)))

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
