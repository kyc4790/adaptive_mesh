function checkEnvelopeTheorem(meshEdge, qFine, dLambda)
    g = getG(meshEdge, qFine);
    disp(dot(g, dLambda));
    disp(dot(g, dLambda) / norm(g) / norm(dLambda));
end