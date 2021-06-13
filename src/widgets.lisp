(in-package :deliscope)

(defclass parsers ()
  ((parsers :initform (make-hash-table :test 'eq :thread-safe t)
            :initarg :parsers
            :accessor parsers))
  )


(defun count-sequences-in-file (file-name &key num)
  (unless (probe-file file-name)
    (error "Could not open file: ~a~%" file-name))
  (let* ((seq-file (sa:make-seq-file-in file-name))
         (title (sa:make-string :char-string))
         (seq (sa:make-string :dna5q-string)))
    (let ((sequences (loop named read-loop
                           for index from 0
                           for res = (progn
                                   (when (sa:at-end seq-file) (return-from read-loop index))
                                   (sa:read-record title seq seq-file))
                           unless (and num (< num 1000000))
                             do (when (and (> index 0) (= 0 (mod index 1000000)))
                                  (format t "~CRead ~a ~a" #\return index file-name)
                                  (finish-output))
                           when (and num (>= index num))
                             do (return-from read-loop index)
                           finally (return-from read-loop index))))
    (sa:close seq-file)
    sequences)))

(defun create-sequence-parser (name files &key num pass-file-name fail-file-name output-file-name overwrite)
  (unless pass-file-name (setf pass-file-name (make-pathname :name (format nil "~a-pass" (string-downcase (string name)))
                                                   :type "fastq")))
  (unless fail-file-name (setf fail-file-name (make-pathname :name (format nil "~a-fail" (string-downcase (string name)))
                                                             :type "fastq")))
  (unless output-file-name (setf output-file-name (make-pathname :name (format nil "~a-results" (string-downcase (string name)))
                                                                 :type "dat")))
  (unless overwrite
    (when (probe-file output-file-name)
      (error "Saving results will fail because the file ~a already exists" output-file-name)))
  (let* ((num-seqs (loop for file in files
                         for num-sequences = (count-sequences-in-file file :num num)
                         do (format t "~CTotal ~a sequences in ~a~%" #\return num-sequences file)
                         collect num-sequences))
         (parser (make-instance 'parser :name name
                                        :files files
                                        :num-sequences-per-file num-seqs
                                        :pass-file-name pass-file-name
                                        :fail-file-name fail-file-name
                                        :output-file-name output-file-name
                                        :overwrite overwrite)))
    (finish-output)
    parser))


(defparameter *progress* nil)
(defparameter *parser* nil)

(defclass progress ()
  ((progress-bar :initarg :progress-bar :reader progress-bar)
   (messages :initarg :messages :reader messages)
   (holder :initarg :holder :accessor holder)))

(defun parse (parser)
  (let* ((stop-button (jw:make-button :description "Stop"))
         (progress-bar (jw:make-int-progress))
         (messages (jw:make-text-area :rows 2 :layout (jw:make-layout :width "300pt"
                                                                      :height "70pt")))
         (holder (jw:make-v-box :children (list progress-bar
                                                messages
                                                stop-button)))
         (progress (make-instance 'progress :progress-bar progress-bar
                                            :messages messages
                                            :holder holder)))
    (let ((thread (bordeaux-threads:make-thread
                   (lambda ()
                     (serial-analyze *parser*
                                     :progress-callback (lambda (progress msg &key wrote-results-to)
                                                          (if wrote-results-to
                                                              (setf (jw:widget-value (messages *progress*))
                                                                    (format nil "~aWrote results to: ~a~%" (jw:widget-value (messages *progress*))
                                                                            wrote-results-to))
                                                              (progn
                                                                (setf (jw:widget-value (progress-bar *progress*)) progress)
                                                                (setf (jw:widget-value (messages *progress*)) msg))))))
                   :name (format nil "parse-~a" (name parser))
                   :initial-bindings (list (cons '*progress* progress)
                                           (cons '*parser* parser)))))
      (jw:on-button-click stop-button (lambda (&rest args)
                                        (declare (ignore args))
                                        (bordeaux-threads:destroy-thread thread)
                                        (setf (jw:widget-value messages) "Stopped parsing.")
                                        )))
    (holder progress)))


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
