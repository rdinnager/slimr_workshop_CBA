rotate.matrix_2 <- function (x, angle = 10, method = "bilinear") 
{
    # if (length(dim(x)) > 2) {
    #     warning("data must be grayscale image")
    # }
    # else {
        # if (max(x) <= 1) {
        #     x <- x * (2^attr(x, "bits.per.sample") - 1)
        # }
        img <- x
        angle.rad <- angle * pi/180
        co.x <- matrix(rep(-(ncol(img)/2 - 0.5):(ncol(img)/2 - 
            0.5), nrow(img)), nrow = nrow(img), byrow = T)
        co.y <- matrix(rep(-(nrow(img)/2 - 0.5):(nrow(img)/2 - 
            0.5), ncol(img)), ncol = ncol(img))
        # if (method == "simple") {
            co.xn <- round(co.x * cos(angle.rad) - co.y * sin(angle.rad))
            co.yn <- round(co.x * sin(angle.rad) + co.y * cos(angle.rad))
            co.xn2 <- co.xn + max(co.xn) + 1
            co.yn2 <- co.yn + max(co.yn) + 1
            img.rot <- numeric(max(co.yn2) * max(co.xn2))
            img.rot[(co.xn2 - 1) * max(co.yn2) + co.yn2] <- img
        # }
        # else if (method == "NN" || method == "bilinear") {
        #     if (ncol(img)%%2 == 0) {
        #         co.xn <- trunc(co.x * cos(angle.rad) - co.y * 
        #           sin(angle.rad) + 0.5)
        #     }
        #     else {
        #         co.xn <- trunc(co.x * cos(angle.rad) - co.y * 
        #           sin(angle.rad))
        #     }
        #     if (nrow(img)%%2 == 0) {
        #         co.yn <- trunc(co.x * sin(angle.rad) + co.y * 
        #           cos(angle.rad) + 0.5)
        #     }
        #     else {
        #         co.yn <- trunc(co.x * sin(angle.rad) + co.y * 
        #           cos(angle.rad))
        #     }
        #     co.list1 <- c(0, 0)
        #     for (i in 1:max(co.xn)) {
        #         y1 <- co.yn[which(co.xn == i)]
        #         y2 <- min(y1):max(y1)
        #         y3 <- matrix(c(rep(i, length(y2)), y2), nrow = 2, 
        #           byrow = T)
        #         co.list1 <- cbind(co.list1, y3)
        #     }
        #     co.list1 <- co.list1[, -1]
        #     co.list2 <- matrix(0, nrow(co.list1), ncol(co.list1))
        #     co.list2[1, ] <- -co.list1[1, ]
        #     if (ncol(img)%%2 == 0) {
        #         co.list2[1, ] <- co.list2[1, ] + 1
        #     }
        #     co.list2[2, ] <- -co.list1[2, ]
        #     if (nrow(img)%%2 == 0) {
        #         co.list2[2, ] <- co.list2[2, ] + 1
        #     }
        #     co.list <- cbind(co.list1, co.list2)
        #     if (ncol(img)%%2 != 0) {
        #         y0 <- co.yn[which(co.xn == 0)]
        #         y02 <- min(y0):max(y0)
        #         y03 <- matrix(c(rep(0, length(y02)), y02), nrow = 2, 
        #           byrow = T)
        #         co.list <- cbind(co.list, y03)
        #     }
        #     co.xn2 <- co.list[1, ] + max(co.list[1, ]) + 1
        #     co.yn2 <- co.list[2, ] + max(co.list[2, ]) + 1
        #     if (ncol(img)%%2 == 0) {
        #         co.list[1, ] <- co.list[1, ] - 0.5
        #     }
        #     if (nrow(img)%%2 == 0) {
        #         co.list[2, ] <- co.list[2, ] - 0.5
        #     }
        #     if (method == "NN") {
        #         if (ncol(img)%%2 == 0) {
        #           co.xn.b <- round(co.list[1, ] * cos(-angle.rad) - 
        #             co.list[2, ] * sin(-angle.rad) + 0.5) + ncol(img)/2
        #         }
        #         else {
        #           co.xn.b <- round(co.list[1, ] * cos(-angle.rad) - 
        #             co.list[2, ] * sin(-angle.rad)) + ncol(img)/2 + 
        #             0.5
        #         }
        #         if (nrow(img)%%2 == 0) {
        #           co.yn.b <- round(co.list[1, ] * sin(-angle.rad) + 
        #             co.list[2, ] * cos(-angle.rad) + 0.5) + nrow(img)/2
        #         }
        #         else {
        #           co.yn.b <- round(co.list[1, ] * sin(-angle.rad) + 
        #             co.list[2, ] * cos(-angle.rad)) + nrow(img)/2 + 
        #             0.5
        #         }
        #         img.rot <- numeric(max(co.yn2) * max(co.xn2))
        #         img.rot[(co.xn2 - 1) * max(co.yn2) + co.yn2] <- img[(co.xn.b - 
        #           1) * max(co.yn.b) + co.yn.b]
        #     }
        #     else if (method == "bilinear") {
        #         co.xn.b <- co.list[1, ] * cos(-angle.rad) - co.list[2, 
        #           ] * sin(-angle.rad) + 0.5 + ncol(img)/2
        #         co.yn.b <- co.list[1, ] * sin(-angle.rad) + co.list[2, 
        #           ] * cos(-angle.rad) + 0.5 + nrow(img)/2
        #         co.x0 <- floor(co.xn.b) + 1
        #         co.y0 <- floor(co.yn.b) + 1
        #         ax <- co.xn.b - co.x0 + 1
        #         ay <- co.yn.b - co.y0 + 1
        #         img.e <- cbind(img[, 1], img, img[, ncol(img)])
        #         img.e <- rbind(img.e[1, ], img.e, img.e[nrow(img.e), 
        #           ])
        #         img.rot <- numeric(max(co.yn2) * max(co.xn2))
        #         img.rot[(co.xn2 - 1) * max(co.yn2) + co.yn2] <- (1 - 
        #           ax) * (1 - ay) * img.e[(co.x0 - 1) * nrow(img.e) + 
        #           co.y0] + ax * (1 - ay) * img.e[co.x0 * nrow(img.e) + 
        #           co.y0] + ay * (1 - ax) * img.e[(co.x0 - 1) * 
        #           nrow(img.e) + co.y0 + 1] + ax * ay * img.e[co.x0 * 
        #           nrow(img.e) + co.y0 + 1]
        #     }
        # }
        dim(img.rot) <- c(max(co.yn2), max(co.xn2))
        attr(img.rot, "bits.per.sample") <- attr(img, "bits.per.sample")
        attr(img.rot, "samples.per.pixel") <- attr(img, "samples.per.pixel")
        return(img.rot)
    # }
}
