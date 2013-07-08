(defpackage :complex-plane-geometry
  (:nicknames :cgeom)
  (:shadow :transformation)
  (:shadowing-import-from :moebius-transformations #:line)
  (:use :clim :clim-lisp :ol :iterate
        :moebius-transformations)
  (:export))

(in-package :complex-plane-geometry)

(define-application-frame complex-plane ()
  ((bound :initform 1000
          :accessor bound)
   (zoom :initform 1
         :accessor zoom)
   (objects :initarg :objects
            :initform nil
            :accessor objects))
  (:panes (canvas :application
                  :width 1000 :height 800
                  :incremental-redisplay t
                  :scroll-bars t)
          (int :interactor
               :width 400 :height 800))
  (:layouts (default (horizontally () int canvas))))

(defparameter transformation
  (make-3-point-transformation*
   ;; from
   0 0
   1 0
   0 1
   ;; to
   500 400
   510 400
   500 390))

;; this gives a drawing region of [-50,50] x [-40, 40]

(defun complex-plane-interactor ()
  (run-frame-top-level (make-instance 'complex-plane)))

(defmacro with-canvas (&body body)
  `(let ((canvas (get-frame-pane *application-frame* 'canvas)))
     ,@body))

(define-presentation-type projective-line ())

(define-presentation-type line () :inherit-from 'projective-line)

(define-presentation-type circle () :inherit-from 'projective-line)

(defmethod add-object ((line line))
  ;; determine the edge points
  (mvbind (min max) (intersect-rectangle line :rmin -50 :rmax 50 :imin -40 :imax 40)
    (let ((point-1 (+ (basepoint line) (* min (direction line))))
          (point-2 (+ (basepoint line) (* max (direction line)))))
      (with-canvas
        (with-output-as-presentation (canvas line 'line)
          (draw-line* canvas
                      (realpart point-1) (imagpart point-1)
                      (realpart point-2) (imagpart point-2)
                      :transformation transformation))))))

(defmethod add-object :after (object)
  (pushnew object (objects *application-frame*)))


(defmethod add-object ((circle circle))
  (with-canvas
    (with-output-as-presentation (canvas circle 'circle))
    (draw-circle* (get-frame-pane *application-frame* 'canvas)
                  (realpart (center circle))
                  (imagpart (center circle))
                  (radius circle)
                  :filled nil
                  :transformation transformation)))

(define-complex-plane-command (com-experiment :name "Experiment" :menu t)
    ()
  (draw-line* *standard-output* 0 0 20 20 :transformation transformation)
  (draw-circle* *standard-output* 0 0 10 :filled nil :transformation transformation))

(define-complex-plane-command (com-hover :name "Hover test" :menu t)
    ((projective-line 'projective-line)
     (line 'line)
     (circle 'circle)))

(define-complex-plane-command (com-add-circle
                               :name "Add circle (center & radius)" :menu t)
    ((center 'number :prompt "Center")
     (radius 'real :prompt "Radius"))
  (add-object (make-instance 'circle :center center :radius radius)))

(define-complex-plane-command (com-add-line
                               :name "Add line (point & direction)" :menu t)
    ((basepoint 'number :prompt "Point")
     (direction 'number :prompt "Direction"))
  (add-object (make-instance 'line :basepoint basepoint :direction direction)))


;;; TODO get presentations to work

;;; TODO filled circles,
;;; TODO filled half-planes,
;;; TODO filled triangles
