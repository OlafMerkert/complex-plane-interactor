(defpackage :moebius-transformations
  (:nicknames :sl2z)
  (:shadow :compose)
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
   #:intersect-rectangle
   #:line-segment
   #:circle-segment
   #:start
   #:end
   #:triangle
   #:triangle-edges
   #:fundamental-domain-base))

(in-package :moebius-transformations)

;;; CAREFUL, we use I as constant here (for convenience)
(defconstant i (complex 0 1))

;;; distance functions
(defun dist (point-1 point-2)
  (abs (- point-1 point-2)))

(defun dist^2 (point-1 point-2)
  (abs^2 (- point-1 point-2)))

(defun abs^2 (point)
  (+ (^ (realpart point) 2)
     (^ (imagpart point) 2)))


;;; **********************************************************************
;;; Moebius transformations

(defclass/f moebius-transformation ()
  (a b c d)
  (:documentation "z |-> (a z + b)/(c z + d)"))

(defun det (mt)
  (with-slots (a b c d) mt
    (- (* a d) (* b c))))

(defun compose (mt-1 mt-2)
  "Obtain composite of moebius transformations by matrix
multiplication."
  (with-slots (a b c d) mt-1
    (with-slots ((u a) (v b) (w c) (x d)) mt-2
      (macrolet ((weave (r1 r2 s1 s2) `(+ (* ,r1 ,s1) (* ,r2 ,s2))))
        (mt (weave a b u w) (weave a b v x)
            (weave c d u w) (weave c d v x))))))

(defun inverse (mt)
  (with-slots (a b c d) mt
    (mt d (- b) (- c) a)))

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

(defun moebius-from-points (p1 p2 p3 &optional (q1 0 q1-p) (q2 1 q2-p) (q3 :infinity q3-p))
  "A moebius transformation is specified by the images q of three
points p (fahnentransitiv)."
  (cond ((or q1-p q2-p q3-p)
         (compose (inverse (moebius-from-points q1 q2 q3))
                  (moebius-from-points p1 p2 p3)))
        ((not (distinct-p p1 p2 p3)) (error "can't determine transformation when two points are the same."))
        ;; formula: z-p1 / z-p3  p2-p3 / p2-p1
        ;; it means we send p1 -> 0, p2 -> 1 and p3 -> oo
        ((eq p1 :infinity)
         (mt (- p2 p3) 0 1 (- p3)))
        ((eq p2 :infinity)
         (mt 1 (- p1) 1 (- p3)))
        ((eq p3 :infinity)
         (mt 1 (- p1) 0 (- p2 p1)))
        (t (let ((f1 (- p2 p3))
                 (f2 (- p2 p1)))
             (mt f1 (* f1 -1 p1)
                 f2 (* f2 -1 p3))))))

;;; **********************************************************************
;;; some complex plane geometric objects
(defclass projective-line ()
  ())

(defclass/f segment ()
  (start
   end))

(defgeneric transform (transformation object))

(defgeneric intersect (object-1 object-2))

;;; **********************************************************************
;;; lines
(defclass/f line (projective-line)
  (basepoint
   direction))

(defclass line-segment (line segment)
  ())

(defun line (basepoint direction)
  (make-instance 'line :basepoint basepoint :direction direction))

(defmethod print-object ((line line) stream)
  (format stream "(line ~A ~A)" (basepoint line) (direction line)))

(defgeneric normalise-direction (line)
  (:documentation "Normalise the `direction' vector to length 1."))

(defmethod normalise-direction ((line line))
  (with-slots (direction) line
    (setf direction
          (/ direction (abs direction)))))

(defmethod normalise-direction ((line-segment line-segment))
  (with-slots (direction start end) line-segment
    (let ((scale (abs direction)))
      (unless (= scale 1)
        (setf direction (/ direction scale)
              start (* start scale)
              end (* end scale))))))

(defmethod initialize-instance :after ((line line) &key)
  (normalise-direction line))

;;; **********************************************************************
;;; circles
(defclass/f circle (projective-line)
  (center
   radius))

(defclass/f circle-segment (circle)
  ;; 1 corresponds to full circle, going counter-clockwise
  (start end))

(defun circle (center radius)
  (make-instance 'circle :center center :radius radius))

(defmethod print-object ((circle circle) stream)
  (format stream "(circle ~A ~A)" (center circle) (radius circle)))

;;; **********************************************************************
;;; triangles

(defclass/f triangle ()
  (vertices)
  (:documentation "represent a triangle with 0 angles at the
  vertices."))

(defparameter fundamental-domain-base
  (make-instance 'circle-segment
                 :center 1/2
                 :radius 1/2
                 :start 0
                 :end 1/2)
  "The base of the 'triangle' with vertices 0,1 and oo.")

(defun line-segment-between-points (point-1 point-2 point-3)
  ""
  ;; find the transformation that sends 0 -> point-1, 1 -> point-2 and
  ;; oo -> point-3
  (let ((mt (inverse (moebius-from-points point-1 point-2 point-3))))
    (transform mt fundamental-domain-base)))

(defun triangle-edges (triangle)
  (dbind (p1 p2 p3) (vertices triangle)
         (list (line-segment-between-points p1 p2 p3)
               (line-segment-between-points p2 p3 p1)
               (line-segment-between-points p3 p1 p2))))



(defun mittelsenkrechte (point-1 point-2)
  (when (equal point-1 point-2)
    (error "no unique line through one point."))
  (make-instance 'line
                 :basepoint (/ (+ point-1 point-2) 2)
                 :direction (* i (- point-1 point-2))))

;;; **********************************************************************
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
        ;; look at the solutions of the quadratic equation.
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

(defmethod translate ((offset number) (line-segment line-segment))
  (make-instance 'line-segment
                 :basepoint (+ offset (basepoint line-segment))
                 :direction (direction line-segment)
                 :start (start line-segment)
                 :end (end line-segment)))

(defmethod translate ((offset number) (circle circle))
  (make-instance 'circle
                 :center (+ offset (center circle))
                 :radius (radius circle)))

(defmethod translate ((offset number) (circle-segment circle-segment))
  (make-instance 'circle-segment
                 :center (+ offset (center circle-segment))
                 :radius (radius circle-segment)
                 :start (start circle-segment)
                 :end (end circle-segment)))

;;; for scaling, we assume that `factor' is never 0.
(defmethod scale ((factor number) (line line))
  (make-instance 'line
                 :basepoint (* factor (basepoint line))
                 :direction (* factor (direction line))))

(defmethod scale ((factor number) (line-segment line-segment))
  (make-instance 'line-segment
                 :basepoint (* factor (basepoint line-segment))
                 :direction (* factor (direction line-segment))
                 :start (start line-segment)
                 :end (end line-segment)))

(defun m+ (&rest args)
  (mod (apply #'+ args) 1))

(defun angle (number)
  (/ (atan number) 2 pi))

(defmethod scale ((factor number) (circle circle))
  (make-instance 'circle
                 :center (* factor (center circle))
                 :radius (* (abs factor) (radius circle))))

(defmethod scale ((factor number) (circle-segment circle-segment))
  (let ((angle (angle factor)))
    (make-instance 'circle-segment
                   :center (* factor (center circle-segment))
                   :radius (* (abs factor) (radius circle-segment))
                   :start (m+ (start circle-segment) angle)
                   :end (m+ (end circle-segment) angle ))))

;;; the (geometrically) most interesting part is inversion at the circle
(defmethod circle-inversion ((circle circle))
  (cond
    ;; if centered at origin, only the radius changes
    ((zerop (center circle))
     (make-instance 'circle :center 0 :radius (/ (radius circle))))
    ;; if origin is on the circle, we get a line
    ((element-p 0 circle)
     (make-instance 'line
                    :basepoint (/ (* 2 (center circle)))
                    :direction (* i (center circle))))
    ;; otherwise we get a circle
    (t (mvbind (point-1 point-2) (intersect (make-instance 'line :basepoint 0 :direction (center circle))
                                            circle)
         ;; we compute the image of the two points on a line through
         ;; the origin, which contains the center
         (let* ((point-1 (/ point-1)) (point-2 (/ point-2))
                (center (/ (+ point-1 point-2) 2)))
           (make-instance 'circle :center center
                          :radius (dist center point-1)))))))

(defmethod circle-inversion ((line line))
  (cond
    ;; if going through the origin, line is just reflected along real axis
    ((element-p 0 line)
     (make-instance 'line :basepoint 0 :direction (/ (direction line))))
    ;; otherwise, we get a circle through the origin
    (t
     (let* ((lotpunkt (intersect line (make-instance 'line :basepoint 0
                                                     :direction (* i (direction line)))))
            (center (/ (* lotpunkt 2))))
      (make-instance 'circle :center center :radius (abs center))))))

;;; TODO try to obtain nice formulas for circles and lines (if
;;; possible)
;;; TODO transformations of segments of circles and lines

;;; **********************************************************************
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
