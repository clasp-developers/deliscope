(in-package :deliscope)


(defun my-as-matrix (sheet &optional na-string)
  "Creates an array from a list of cells of the form ((:A . 1) . 42)
Empty columns or rows are ignored (column and row names are returned as additional values).
When a value is equal to na-string, nil is returned instead."
  (declare (optimize debug ))
  (let* ((refs (mapcar #'first sheet))
         (cols (sort (remove-duplicates (mapcar #'car refs)) #'string< :key (lambda (x) (format nil "~3@A" x))))
         (rows (sort (remove-duplicates (mapcar #'cdr refs)) #'<))
         (last-row (car (last rows)))
         (output (make-array (list last-row (length cols)) :initial-element nil)))
    (loop for ((col . row) . val) in sheet
       do (setf (aref output (1- row) (the fixnum (position col cols)) )
                (if (equal val na-string) nil val)))
    (values output cols rows)))


(defun parse-plate (map result tlx tly &optional (brx (+ tlx 12)) (bry (+ tly 8))  (offset 13))
  (loop for x from tlx below brx
        do (loop for y from tly below bry
                 for code = (aref map y x)
                 for val = (aref map y (- x offset))
                 for code-no-slash = (format nil "~a~a" (subseq code 0 4) (subseq code 5 9))
                 do (setf (gethash code result) val)
                 do (setf (gethash code-no-slash result) val)
                 do (format t "~a = ~a~%" code val))))

(defun parse-map (m)
  (let ((ht (make-hash-table :test 'equalp)))
    (parse-plate m ht 14 2)
    (parse-plate m ht 14 15)
    (parse-plate m ht 14 28)
    ht))
    
(defun read-map (filename)
  (let* ((sheet (xlsx:read-sheet filename 1))
         (data (my-as-matrix sheet))
         (dims (array-dimensions data))
         (rows (elt dims 0))
         (cols (elt dims 1)))
    data))

(defun translate (m seq)
  (loop for mon in seq
        for val = (gethash mon m)
        do (format t "~a -> ~a~%" mon val)
        collect val))
