(defpackage :complex-plane-geometry
  (:nicknames :cgeom)
  (:shadow :transformation)
  (:shadowing-import-from :moebius-transformations #:line)
  (:use :clim :clim-lisp :ol :iterate
        :moebius-transformations)
  (:export))

(in-package :complex-plane-geometry)

(defparameter window-width 1000)
(defparameter window-height 800)

(define-application-frame complex-plane ()
  ((bound :initform 1000
          :accessor bound)
   (zoom :initform 1
         :accessor zoom)
   (objects :initarg :objects
            :initform nil
            :accessor objects))
  (:panes (canvas :application
                  :width window-width :height window-height
                  :incremental-redisplay t
                  :scroll-bars t)
          (int :interactor
               :width 400 :height window-height))
  (:layouts (default (horizontally () int canvas))))

(progn
  (defparameter pixels-per-one 100)

  (let ((h-center (floor window-width 2))
        (v-center (floor  window-height 2)))
    (defparameter transformation
      (make-3-point-transformation*
       ;; from
       0 0
       1 0
       0 1
       ;; to
       h-center v-center
       (+ h-center pixels-per-one) v-center
       h-center (- v-center pixels-per-one)))))

;; this gives a drawing region of [-50,50] x [-40, 40]

(defun complex-plane-interactor ()
  (run-frame-top-level (make-instance 'complex-plane)))

(defmacro with-canvas (&body body)
  `(let ((canvas (get-frame-pane *application-frame* 'canvas)))
     ,@body))

(define-presentation-type projective-line ())

(define-presentation-type line () :inherit-from 'projective-line)
(define-presentation-type line-segment () :inherit-from 'line)

(define-presentation-type circle () :inherit-from 'projective-line)
(define-presentation-type circle-segment () :inherit-from 'circle)

(define-presentation-type triangle ())

(defmethod add-object ((line line))
  ;; determine the edge points
  (mvbind (min max) (intersect-rectangle line :rmin -50 :rmax 50 :imin -40 :imax 40)
    (let ((point-1 (+ (basepoint line) (* min (direction line))))
          (point-2 (+ (basepoint line) (* max (direction line)))))
      (with-canvas
        (with-output-as-presentation (canvas line 'line
                                             :allow-sensitive-inferiors nil
                                             :single-box t)
          (draw-line* canvas
                      (realpart point-1) (imagpart point-1)
                      (realpart point-2) (imagpart point-2)
                      :transformation transformation))))))

(defmethod add-object ((line-segment line-segment))
  (with-slots ((min start) (max end)) line-segment
    (let ((point-1 (+ (basepoint line-segment) (* min (direction line-segment))))
          (point-2 (+ (basepoint line-segment) (* max (direction line-segment)))))
      (with-canvas
        (with-output-as-presentation (canvas line-segment 'line-segment
                                             :allow-sensitive-inferiors nil
                                             :single-box t)
          (draw-line* canvas
                      (realpart point-1) (imagpart point-1)
                      (realpart point-2) (imagpart point-2)
                      :transformation transformation))))))


(defmethod add-object :after (object)
  (pushnew object (objects *application-frame*)))


(defmethod add-object ((circle circle))
  (with-canvas
    (with-output-as-presentation (canvas circle 'circle
                                         :allow-sensitive-inferiors nil
                                         :single-box t)
      (draw-circle* canvas
                    (realpart (center circle))
                    (imagpart (center circle))
                    (radius circle)
                    :filled nil
                    :transformation transformation))))

(defmethod add-object ((circle-segment circle-segment))
  (with-canvas
    (with-output-as-presentation (canvas circle-segment 'circle-segment
                                         :allow-sensitive-inferiors nil
                                         :single-box t)
      (draw-circle* canvas
                    (realpart (center circle-segment))
                    (imagpart (center circle-segment))
                    (radius circle-segment)
                    :start-angle (* 2 pi (end circle-segment))
                    :end-angle   (* 2 pi (start circle-segment))
                    :filled nil
                    :transformation transformation))))

(defmethod add-object ((triangle triangle))
  (with-canvas
    (with-output-as-presentation (canvas triangle 'triangle
                                         :allow-sensitive-inferiors t
                                         :single-box t)
      (mapc #'add-object (triangle-edges triangle)))))

(define-complex-plane-command (com-quit :name "Quit") ()
  (frame-exit *application-frame*))

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

(define-complex-plane-command (com-add-triangle
                               :name "Add triangle (vertices)" :menu t)
    ((point-1 'number :prompt "Vertex 1")
     (point-2 'number :prompt "Vertex 2")
     (point-3 'number :prompt "Vertex 3"))
  (add-object (make-instance 'triangle :vertices (list point-1 point-2 point-3))))

(define-complex-plane-command (com-add-3points
                               :name "Add projective line through 3 points" :menu t)
    ((point-1 'number :prompt "Point 1")
     (point-2 'number :prompt "Point 2")
     (point-3 'number :prompt "Point 3"))
    (add-object (projective-line-3points point-1 point-2 point-3)))

(define-complex-plane-command (com-translate
                               :name "Translate" :menu t)
    ((object 'projective-line) (offset 'number))
  (add-object (translate offset object)))

(define-complex-plane-command (com-scale
                               :name "Scale" :menu t)
    ((object 'projective-line) (factor 'number))
  (add-object (scale factor object)))

(define-complex-plane-command (com-invert
                               :name "Invert" :menu t)
    ((object 'projective-line))
  (add-object (circle-inversion object)))

(define-complex-plane-command (com-mittelsenkrechte
                               :name t :menu nil)
    ((point-1 'number :prompt "Point 1")
     (point-2 'number :prompt "Point 2"))
  (add-object (mittelsenkrechte point-1 point-2)))

(define-complex-plane-command (com-draw-farey-sets :name "Draw Farey sets" :menu t)
    ()
  (dolist (r ccfg:regions)
    (unless (symbolp r)
      (add-object  r))))

(define-complex-plane-command (com-clear-canvas :name "Clear Canvas" :menu t)
    ()
  (with-canvas (window-clear canvas)))

(define-complex-plane-command (com-test-stuff :name "Test stuff" :menu t)
    ()
  #|(add-object fundamental-domain-base)|#
  #|(add-object (make-instance 'line-segment :basepoint 0 :direction (+ 1 i)
  :start -1 :end 2))|#
  (add-object (make-instance 'triangle :vertices (list 0 1 i))))

;;; TODO get presentations to work

;;; TODO filled circles,
;;; TODO filled half-planes,
;;; TODO filled triangles
