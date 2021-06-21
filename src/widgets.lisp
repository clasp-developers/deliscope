(in-package :deliscope)

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


(defparameter *number-of-threads* (core:num-logical-processors))

(defun set-number-of-threads (num)
  (setf num (max 1 num))
  (setf *number-of-threads* (min num (core:num-logical-processors))))

(defparameter *progress* nil)
(defparameter *parser* nil)

(defclass progress ()
  ((progress-bar :initarg :progress-bar :reader progress-bar)
   (messages :initarg :messages :reader messages)
   (holder :initarg :holder :accessor holder)))

(defun parse (&rest args)
  "Usage: (parse [:limit <integer>] <parsers>...)

Parse the sequences specified by each parser.
You can limit the number sequences loaded from each file within each parser by adding :limit <num>.
"
  (let ((limit nil))
    (when (eq (car args) :limit)
      (pop args)
      (setf limit (pop args))
      (unless (integerp limit)
        (error "limit must be an integer")))
    (let* ((parsers args)
           (container (make-instance 'jw:v-box))
           (panel (cw:make-threaded-task-page
                   container
                   (format nil "Parse")
                   (lambda (instance action parsers progress-callback)
                     (declare (ignore action instance))
                     (format t "Number of threads: ~a~%" *number-of-threads*)
                     (analyze-parsers-using-workers
                      *number-of-threads*
                      parsers
                      :max-sequences-per-file limit
                      :progress-callback
                      (let ((last-val -1))
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
                                    (funcall progress-callback :value val
                                                               :maximum 100))
                                  )
                                )))))
                     t)
                   :parameter parsers
                   :label "Click button to start.")))
      (j:display container)
      (cw:run-task panel)
      (values))))


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
