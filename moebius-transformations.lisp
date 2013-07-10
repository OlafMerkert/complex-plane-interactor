(defpackage :moebius-transformations
  (:nicknames :sl2z)
  (:use :cl :ol :iterate )
  (:export
   #:i
   #:moebius-transformation
   #:evaluate
   #:mt
   #:well-defined-p
   #:projective-line
   #:line
   #:basepoint
   #:direction
   #:circle
   #:center
   #:radius
   #:dist
   #:dist^2
   #:abs^2
   #:mittelsenkrechte
   #:intersect
   #:element-p
   #:distinct-p
   #:collinear-p
   #:projective-line-3points
   #:transform
   #:translate
   #:scale
   #:circle-inversion
   #:det
   #:intersect-rectangle))

(in-package :moebius-transformations)

;;; CAREFUL, we use I as constant here (for convenience)
(defconstant i (complex 0 1))

(defclass/f moebius-transformation ()
  (a b c d)
  (:documentation "z |-> (a z + b)/(c z + d)"))

(defun det (mt)
  (with-slots (a b c d) mt
    (- (* a d) (* b c))))

(defmethod evaluate ((mt moebius-transformation) (z number))
  (with-slots (a b c d) mt
    (/ (+ (* a z) b)
       (+ (* c z) d))))

(defun mt (a b c d)
  (make-instance 'moebius-transformation
                 :a a :b b :c c :d d))

(defun well-defined-p (mt)
  (with-slots (a b c d) mt
    (not (zerop (- (* a d) (* b c))))))

;;; some complex plane geometric objects
(defclass projective-line ()
  ())

(defclass/f line (projective-line)
  (basepoint
   direction))

(defun line (basepoint direction)
  (make-instance 'line :basepoint basepoint :direction direction))

(defmethod print-object ((line line) stream)
  (format stream "(line ~A ~A)" (basepoint line) (direction line)))

(defmethod initialize-instance :after ((line line) &key)
  ;; normalise direction to length 1
  (with-slots (direction) line
    (setf direction
          (/ direction (abs direction)))))

(defclass/f circle (projective-line)
  (center
   radius))

(defun circle (center radius)
  (make-instance 'circle :center center :radius radius))

(defmethod print-object ((circle circle) stream)
  (format stream "(circle ~A ~A)" (center circle) (radius circle)))

(defun dist (point-1 point-2)
  (abs (- point-1 point-2)))

(defun dist^2 (point-1 point-2)
  (abs^2 (- point-1 point-2)))

(defun abs^2 (point)
  (+ (expt (realpart point) 2)
     (expt (imagpart point) 2)))


(defun mittelsenkrechte (point-1 point-2)
  (when (equal point-1 point-2)
    (error "no unique line through one point."))
  (make-instance 'line
                 :basepoint (/ (+ point-1 point-2) 2)
                 :direction (* i (- point-1 point-2))))

;;; some functions for incidence of geometric objects
(defmethod intersect ((line-1 line) (line-2 line))
  (let ((a (realpart (direction line-1)))
        (c (imagpart (direction line-1)))
        (b (realpart (direction line-2)))
        (d (imagpart (direction line-2))))
    (let ((det (- (* a d) (* b c))))
      (cond ((not (zerop det))
             ;;  we have a unique intersection
             (let* ((dif (- (basepoint line-1) (basepoint line-2)))
                    (x (realpart dif)) (y (imagpart dif)))
               (complex
                (/ (- (* d x) (* c y)) det)
                (/ (- (* a y) (* b x)) det))))
            ;; same lines
            ((element-p (basepoint line-1) line-2)
             line-1)
            ;; parallel
            (t nil)))))

(defmethod intersect ((line line) (circle circle))
  (let ((p (basepoint line))
        (q (center circle)))
    (let ((a (abs^2 (direction line)))
          (b (* 2 (realpart (* (direction line) (conjugate (- p q))))))
          (c (- (abs^2 (- p q)) (expt (radius circle) 2))))
      (let* ((d (sqrt (- (expt b 2) (* 4 a c))))
             (l1 (/ (+ (- b) d) 2 a))
             (l2 (/ (- (- b) d) 2 a)))
        (cond ((zerop d) (+ p (* l1 (direction line))))
              ((realp d) (values (+ p (* l1 (direction line)))
                                 (+ p (* l2 (direction line)))))
              (t nil))))))

(defmethod intersect ((circle circle) (line line))
  (intersect line circle))

(defmethod intersect ((circle-1 circle) (circle-2 circle))
  (unless (< (+ (radius circle-1) (radius circle-2))
             (dist (center circle-1) (center circle-2)))
    ;; TODO
    ))

(defmethod element-p ((point number) (line line))
  (realp (/ (- point (basepoint line))
            (direction line))))

(defmethod element-p ((point number) (circle circle))
  (= (dist^2 point (center circle))
     (expt (radius circle) 2)))


(defun distinct-p (arg &rest args)
  "Check that all arguments are different."
  (if (null args) t
      (when (every (complement (lambda (x) (equal x arg))) args)
        (apply #'distinct-p args))))

(defun collinear-p (point-1 point-2 point-3)
  (let ((dir2 (- point-2 point-1))
        (dir3 (- point-3 point-1)))
    ;; assume not all points are the same for now
    (cond ((zerop dir2) dir3)
          ((realp (/ dir3 dir2)) dir2)
          (t nil))))

(defun projective-line-3points (point-1 point-2 point-3)
  "3 distinct points on the complex plane either define a line or a circle"
  (unless (and (numberp point-1)
               (numberp point-2)
               (numberp point-3))
    (error "Expected complex numbers"))
  (unless (distinct-p point-1 point-2 point-3)
    (error "Expected distinct points."))
  (aif (collinear-p point-1 point-2 point-3)
       (make-instance 'line :basepoint point-1 :direction it)
       (let ((center (compute-circle-center-3points point-1 point-2 point-3)))
         (make-instance 'circle :center center :radius (dist point-1 center)))))

(defun compute-circle-center-3points (point-1 point-2 point-3)
  (let ((ms2 (mittelsenkrechte point-1 point-2))
        (ms3 (mittelsenkrechte point-1 point-3)))
    (intersect ms2 ms3)))

;; moebius transformations map projective lines on projective lines.
(defmethod transform ((mt moebius-transformation) (line projective-line))
  (with-slots (a b c d) mt
    (if (zerop c)
        ;; TODO direct implementation (perhaps)
        (translate (/ b d) (scale (/ a d) line))
        ;;
        (let* ((det (det mt))
               (a* (/ (- c) det))
               (b* (/ a det)))
          (translate b* (scale a* (circle-inversion (translate (/ d c) line))))))))

(defmethod translate ((offset number) (line line))
  (make-instance 'line
                 :basepoint (+ offset (basepoint line))
                 :direction (direction line)))

(defmethod translate ((offset number) (circle circle))
  (make-instance 'circle
                 :center (+ offset (center circle))
                 :radius (radius circle)))

(defmethod scale ((factor number) (line line))
  (make-instance 'line
                 :basepoint (* factor (basepoint line))
                 :direction (* factor (direction line))))

(defmethod scale ((factor number) (circle circle))
  (make-instance 'circle
                 :center (* factor (center circle))
                 :radius (* (abs factor) (radius circle))))

(defmethod circle-inversion ((circle circle))
  (cond
    ;; if centered at origin, nothing happens
    ((zerop (center circle))
     circle)
    ;; if origin is on the circle, we get a line
    ((element-p 0 circle)
     (make-instance 'line
                    :basepoint (/ (* 2 (center circle)))
                    :direction (* i (center circle))))
    ;; otherwise we get a circle
    (t (mvbind (point-1 point-2) (intersect (make-instance 'line :basepoint 0 :direction (center circle))
                                            circle)
         (let* ((point-1 (/ point-1)) (point-2 (/ point-2))
                (center (/ (+ point-1 point-2) 2)))
           (make-instance 'circle :center center
                          :radius (dist center point-1)))))))

(defmethod circle-inversion ((line line))
  (cond
    ;; if going through the origin, line is just reflected along real axis
    ((element-p 0 line)
     (make-instance 'line :basepoint 0
                    :direction (/ (direction line))))
    ;; otherwise, we get a circle through the origin
    (t
     (let* ((lotpunkt (intersect line (make-instance 'line :basepoint 0
                                                     :direction (* i (direction line)))))
            (center (/ (* lotpunkt 2))))
      (make-instance 'circle :center center :radius (abs center))))))

;;; TODO try to obtain nice formulas for circles and lines (if possible)

;;; intersections of lines with rectangles
(defun intersect-rectangle (line &key rmin rmax imin imax)
  ;; TODO use a sorting network (have always <=4 entries)
  (let ((ls (sort (append (endpoints-line-real line rmin rmax)
                          (endpoints-line-imag line imin imax))
                  #'<=)))
    (values-list (if (< 2 (length ls))
                     (subseq ls 1 3)
                     ls))))

(defun endpoints (p v min max)
  (list (/ (- min p) v)
        (/ (- max p) v)))

(bind-multi ((endpoints-line-real endpoints-line-real endpoints-line-imag)
             (realpart realpart imagpart))
  (defun endpoints-line-real (line min max)
                 (let ((v (realpart (direction line))))
                   (unless (zerop v)
                     (endpoints (realpart (basepoint line)) v min max)))))
