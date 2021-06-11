(in-package :analyze-del)

(defun quality-summary ()
  (loop for qcode in '(#\E #\A #\< #\6 #\/ #\!)
        do (format t "~a -> ~a~%" qcode (- (char-int qcode) (char-int #\!)))))

(defun get-quality-data (seq)
  (loop for pos below (sa:length seq)
        collect (cons pos (sa::get-quality-value (sa::[]& seq pos)))))

(defun filter-seq-file (seq-file &key (num 10) verbose)
  (let ((title (sa:make-string :char-string))
        (seq (sa:make-string :dna5q-string))
        (count 0)
        (count-sequences-without-n 0)
        (low-quality-count 0))
    (when verbose (format t "About to read~%"))
    (loop named read-loop
          for index from 0
          for res = (progn
                      (when (sa:at-end seq-file) (return-from read-loop nil))
                      (incf count)
                      (when (and num (>= count num))
                        (return-from read-loop nil))
                      (when (and (or (null num) (> num 1000000))
                                 (= 0 (mod count 1000000)))
                        (format t "Read ~a sequences and ~a without N  ~a low-quality~%" count count-sequences-without-n low-quality-count))
                      (sa:read-record title seq seq-file)
                      (let ((lqc (sa::count-quality-value-less-than seq 20)))
                        (when (> lqc 1)
                          (incf low-quality-count)))
                      (let ((seq-str (sa:to-string seq)))
                        (unless (position #\N seq-str)
                          (incf count-sequences-without-n)))))
    (format t "Read ~a total sequences and ~a sequences without N ~a low-quality~%"
            count
            count-sequences-without-n low-quality-count)))
