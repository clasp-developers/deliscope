(in-package :deliscope)

(defclass parsers ()
  ((parsers :initform (make-hash-table :test 'eq :thread-safe t)
            :initarg :parsers
            :accessor parsers))
  )

(defun count-sequences-in-file (file-name &key num verbose)
  (declare (ignore verbose))
  (unless (probe-file file-name)
    (error "Could not find file: ~a" file-name))
  (if num
      num
      (let ((num-lines (core:count-lines-in-file file-name 0 2)))
        (/ num-lines 4))))
  
(defun create-sequence-parser (name files &key num pass-file-name fail-file-name output-file-name (overwrite t))
  (unless pass-file-name
    (setf pass-file-name (make-pathname :name (format nil "~a-pass" (string-downcase (string name)))
                                        :type "fastq")))
  (unless fail-file-name
    (setf fail-file-name (make-pathname :name (format nil "~a-fail" (string-downcase (string name)))
                                        :type "fastq")))
  (unless output-file-name
    (setf output-file-name (make-pathname :name (format nil "~a-results" (string-downcase (string name)))
                                          :type "dat")))
  (unless overwrite
    (when (probe-file output-file-name)
      (error "Saving results will fail because the file ~a already exists" output-file-name)))
  (format t "Scanning the files to count the number of sequences.~%")
  (format t "This may take a few minutes - hit ii to interrupt.~%")
  (finish-output)
  (let* ((num-seq-per-file (lparallel:pmapcar (lambda (file)
                                     (count-sequences-in-file file :num num :verbose nil))
                                   files))
         (parser (make-instance 'parser :name name
                                        :files files
                                        :num-sequences-per-file num-seq-per-file
                                        :pass-file-name pass-file-name
                                        :fail-file-name fail-file-name
                                        :output-file-name output-file-name
                                        :overwrite overwrite)))
    (format t "Number of sequences: ~a~%" (apply '+ num-seq-per-file))
    (finish-output)
    parser))


(defparameter *progress* nil)
(defparameter *parser* nil)

(defclass progress ()
  ((progress-bar :initarg :progress-bar :reader progress-bar)
   (messages :initarg :messages :reader messages)
   (holder :initarg :holder :accessor holder)))


(defun run-impl (serial-parallel parsers)
  (let* ((container (make-instance 'jw:v-box))
         (panels (loop for parser in parsers
                       collect (cw:make-threaded-task-page
                                container
                                (format nil "Task ~a" (name parser))
                                (lambda (instance action parser progress-callback)
                                  (declare (ignore action))
                                  (funcall serial-parallel
                                           parser
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
                                                                            (funcall progress-callback :value val :maximum 100))))))))
                                  t)
                                :parameter parser
                                :label "Click button to start."))))
    (j:display container)
    (loop for panel in panels
          do (cw:run-task panel))
    (values)))

(defun run (&rest parsers)
  "Parse sequence files one at a time but parsers are run in parallel."
  (run-impl 'serial-analyze parsers))

(defun run-fast (&rest parsers)
  "Parse sequence files in parallel and parsers in parallel.
 This is experimental and depends on the SeqAn 2.0 library to be thread safe,
 which I am not absolutely certain it is.  If this crashes - use 'run'."
  (run-impl 'parallel-analyze parsers))

(defun plot-counts (data title)
  (let* ((x-data (coerce (subseq (mapcar #'car data) 0) 'vector))
         (y-data (coerce (subseq (mapcar #'cdr data) 0) 'vector)))
    (dv:bar-plot x-data
                 y-data
                 :width 600
                 :x-title "index" :y-title "Log10 counts"
                 :title title)))

(defun plot-unfiltered-counts (analysis)
  (let ((data (sequence-counts (sequences analysis))))
    (plot-counts data (format nil "~a unfiltered counts" (name analysis)))))

(defun plot-filtered-counts (analysis)
  (let ((data (sequence-counts (filtered-sequences analysis))))
    (plot-counts data (format nil "~a filtered counts" (name analysis)))))
